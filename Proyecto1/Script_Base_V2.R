# Cargado de las librerias
rm(list = ls())

library(readxl)
library(Amelia)
library(tidyverse)
library(scales)
library(sf)
library(gstat)
library(ade4)
library(vegan)
library(ggplot2)
library(rnaturalearth)
library(geodata)

# Lectura del dataset


data = read_rds("Data/Temperatura_Ambiente_del_Aire_13_07_2026.rds")

data <- data |> filter(!Departamento %in% c("AMAZONAS","VICHADA"))

# Revision de la base de datos

## Verificar datos missing

missmap(data)

colSums(is.na(data))

## Verificar que las columnas tengan los nombres sin espacios, acentos u otros

names(data)

## Ver los primeros registros

head(data)

## Tamaño de la base de datos

dim(data)

# Agrupado por localización

## Crear la columna de agrupacion

data$Ubicacion <- paste(data$Latitud, "_", data$Longitud)


## Agrupar por la localizacion

data_summary <- data |>
  group_by(Ubicacion) |> 
  mutate(TemperaturaMedia = mean(ValorObservado),
            TemperaturaQ1 = quantile(ValorObservado, 0.25),
            TemperaturaMediana = median(ValorObservado),
            TemperaturaQ3 = quantile(ValorObservado, 0.75),
            Registros = n(),
            Sd = sd(ValorObservado),
            Diff_Median_Mean = TemperaturaMedia - TemperaturaMediana) |> 
  distinct(Ubicacion, .keep_all = TRUE) |> 
  ungroup()

# Filtrado de columnas

names_columnas = c("CodigoEstacion", "ValorObservado", "NombreEstacion", "Departamento", "Municipio",
                   "ZonaHidrografica", "Longitud", "Latitud", "TemperaturaMediana", "Registros")


data_geo_tem = data_summary %>% select(names_columnas)


################################################################################
################## Analisis previos ############################################
################################################################################

summary(data_geo_tem[, c("TemperaturaMediana")])


summary(data_summary[, c("Longitud", "Latitud")])


punto_referencia <- st_sfc(st_point(c(-75, 5)), crs = 4326)

candidatos <- c("MAGNA-SIRGAS / Colombia Bogota zone (3116)" = 3116,
                "MAGNA-SIRGAS 2018 / Origen-Nacional (9377)" = 9377,
                "WGS84 / UTM zona 18N (32618)"                = 32618)


datos_sf <- st_as_sf(data_geo_tem, coords = c("Longitud", "Latitud"), crs = 4326, remove = FALSE) |>
  st_transform(9377)

datos_sf # Longitud y Latitud ya estan en metros



# Mapa por departamento

# Descargar departamentos de Colombia 
depto_colombia <- gadm(country = "COL", level = 1, path = "Data") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

# 2. Graficar mapa por departamentos con puntos en color continuo
ggplot() +
  geom_sf(data = depto_colombia, fill = "white", color = "black") +
  geom_sf(data = datos_sf, aes(color = TemperaturaMediana), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo", name = "Temp. (°C)") +
  coord_sf(crs = 4326) +
  theme_bw() +
  labs(
    title = "Ubicación de Estaciones por Departamento",
    subtitle = "Color continuo: Temperatura (°C)",
    x = "Longitud",
    y = "Latitud"
  )


# Semivariograma cloud

vc_temp <- variogram(TemperaturaMediana ~ 1, datos_sf, cloud = TRUE)

ggplot(as.data.frame(vc_temp), aes(x = dist, y = gamma)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred") +
  theme_test(base_size = 13) +
  labs(x = "Distancia geográfica (m)", y = expression((Z(s[i]) - Z(s[j]))^2),
       title = "Semivariograma cloud — temperatura") +
  theme(axis.title = element_text(face = "bold"))


# Test de mantel

geodist   <- dist(st_coordinates(datos_sf))
dist_temp <- dist(datos_sf$TemperaturaMediana)

set.seed(12345)
mantel_vegan <- vegan::mantel(geodist, dist_temp, permutations = 999)
mantel_vegan

mantel_ade4 <- ade4::mantel.rtest(geodist, dist_temp, nrepet = 999)
mantel_ade4


# Semivariograma empírico

distancia_maxima <-  max(dist(st_coordinates(datos_sf)))/2
ancho_intervalo  <- distancia_maxima / 20

cat("Distancia máxima entre estaciones (completa):",
    round(max(dist(st_coordinates(datos_sf))), 0), "m\n")

cat("Distancia máxima usada (cutoff):", round(distancia_maxima, 0), "m\n")

cat("Ancho de cada intervalo:", round(ancho_intervalo, 0), "m\n")

v_temp <- variogram(TemperaturaMediana ~ 1, datos_sf, cutoff = distancia_maxima, width = ancho_intervalo)
v_temp


ggplot(v_temp, aes(x = dist, y = gamma)) +
  geom_point(aes(size = np), color = "#3a7ca5") +
  geom_line(color = "#3a7ca5", linetype = "dashed") +
  scale_size_continuous(name = "N° pares", range = c(1.5, 5)) +
  theme_test(base_size = 13) +
  labs(x = "Distancia (m)", y = expression(hat(gamma)(h)),
       title = "Semivariograma empírico (sin remover tendencia) — temperatura") +
  theme(axis.title = element_text(face = "bold"))


# Diagonostico de tendencia

data_summary %>%
  pivot_longer(cols = c(Longitud, Latitud), names_to = "coordenada", values_to = "valor") %>%
  ggplot(aes(x = valor, y = TemperaturaMediana)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "darkred") +
  facet_wrap(~coordenada, scales = "free_x") +
  theme_test(base_size = 13) +
  labs(x = NULL, y = "Temperatura (°C)",
       title = "Temperatura vs. coordenadas — evidencia de tendencia") +
  theme(axis.title = element_text(face = "bold"))


# Diagnostico de tendencia removiendo la tendencia en Latitud


data_summary$TemperaturaMedianaReg <- lm(TemperaturaMediana~Latitud, data = data_geo_tem)[["residuals"]]

datos_sf$TemperaturaMedianaReg <- data_summary$TemperaturaMedianaReg

# Test de Mantel removiendo la tendencia

geodist   <- dist(st_coordinates(datos_sf))
dist_temp <- dist(datos_sf$TemperaturaMedianaReg)

set.seed(12345)
mantel_vegan <- vegan::mantel(geodist, dist_temp, permutations = 999)
mantel_vegan

mantel_ade4 <- ade4::mantel.rtest(geodist, dist_temp, nrepet = 999)
mantel_ade4


data_summary %>%
  pivot_longer(cols = c(Longitud, Latitud), names_to = "coordenada", values_to = "valor") %>%
  ggplot(aes(x = valor, y = TemperaturaMedianaReg)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "darkred") +
  facet_wrap(~coordenada, scales = "free_x") +
  theme_test(base_size = 13) +
  labs(x = NULL, y = "Temperatura (°C)",
       title = "Temperatura vs. coordenadas — evidencia de tendencia") +
  theme(axis.title = element_text(face = "bold"))

# Normalidad

# Sin tendencia 

v_1er_orden <- variogram(TemperaturaMediana ~ Latitud, datos_sf,
                         cutoff = distancia_maxima, width = ancho_intervalo)

v_2do_orden <- variogram(TemperaturaMediana ~ poly(Latitud, degree = 2), datos_sf,
                         cutoff = distancia_maxima, width = ancho_intervalo)

bind_rows(v_temp      %>% mutate(modelo = "Sin remover tendencia"),
          v_1er_orden %>% mutate(modelo = "Tendencia 1er orden removida"),
          v_2do_orden %>% mutate(modelo = "Tendencia 2do orden removida")) %>%
  ggplot(aes(x = dist, y = gamma, color = modelo)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("grey40", "#3a7ca5", "#d9534f"), name = NULL) +
  theme_test(base_size = 13) +
  labs(x = "Distancia (m)", y = expression(hat(gamma)(h)),
       title = "Efecto de remover la tendencia sobre el semivariograma") +
  theme(legend.position = "bottom", axis.title = element_text(face = "bold"))

# Semivariogramas direccionales

v_direccional <- variogram(TemperaturaMediana ~ Latitud, datos_sf,
                           alpha = c(0, 45, 90, 135), cutoff = distancia_maxima)

ggplot(v_direccional, aes(x = dist, y = gamma, color = factor(dir.hor))) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(name = "Dirección (°)") +
  theme_test(base_size = 13) +
  labs(x = "Distancia (m)", y = expression(hat(gamma)(h)),
       title = "Semivariogramas direccionales — temperatura (tendencia removida)") +
  theme(axis.title = element_text(face = "bold"))

# Envolvente de aleatorizacion Monte Carlo

set.seed(12345)
n_simulaciones <- 199

envolvente_sim <- matrix(NA_real_, nrow = nrow(v_temp), ncol = n_simulaciones)
datos_sf_perm  <- datos_sf

for (i in seq_len(n_simulaciones)) {
  datos_sf_perm$temp_perm <- sample(datos_sf$TemperaturaMediana)
  v_perm <- variogram(temp_perm ~ 1, datos_sf_perm,
                      cutoff = distancia_maxima, width = ancho_intervalo)
  envolvente_sim[, i] <- v_perm$gamma
}

envolvente <- tibble(dist  = v_temp$dist,
                     gamma = v_temp$gamma,
                     inf   = apply(envolvente_sim, 1, quantile, probs = 0.025),
                     sup   = apply(envolvente_sim, 1, quantile, probs = 0.975)) %>%
  mutate(fuera_envolvente = gamma < inf | gamma > sup)

envolvente

ggplot(envolvente, aes(x = dist)) +
  geom_ribbon(aes(ymin = inf, ymax = sup), fill = "grey80", alpha = 0.6) +
  geom_line(aes(y = gamma), color = "#d9534f", linewidth = 1) +
  geom_point(aes(y = gamma), color = "#d9534f", size = 2) +
  theme_test(base_size = 13) +
  labs(x = "Distancia (m)", y = expression(hat(gamma)(h)),
       title = "Semivariograma observado (rojo) vs. envolvente de aleatorización (gris)") +
  theme(axis.title = element_text(face = "bold"))

mean(envolvente$fuera_envolvente)

# Normalidad de la variable

ggplot(data_geo_tem, aes(x = TemperaturaMediana)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12, fill = "#3a7ca5", color = "white") +
  geom_density(color = "darkred", linewidth = 1) +
  theme_test(base_size = 13) +
  labs(x = "Temperatura (°C)", y = "Densidad", title = "Histograma") +
  theme(axis.title = element_text(face = "bold"))

ggplot(data_geo_tem, aes(sample = TemperaturaMediana)) +
  stat_qq() +
  stat_qq_line(color = "darkred") +
  theme_test(base_size = 13) +
  labs(x = "Cuantiles teóricos (Normal)", y = "Cuantiles muestrales",
       title = "Gráfico cuantil-cuantil (Q-Q)") +
  theme(axis.title = element_text(face = "bold"))

shapiro.test(data_geo_tem$TemperaturaMediana)

asimetria_temp <- mean((data_geo_tem$TemperaturaMediana - mean(data_geo_tem$TemperaturaMediana))^3) /
  sd(data_geo_tem$TemperaturaMediana)^3
asimetria_temp

