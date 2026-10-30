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


# Lectura del dataset


data = read_rds("Data/Temperatura_Ambiente_del_Aire_13_07_2026.rds")



# Revision de la base de datos

## Verificar datos missing

missmap(data)

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


data_geo_tem = data_summary |> select(names_columnas)


################################################################################
################## Analisis previos ############################################
################################################################################

summary(data_geo_tem[, c("TemperaturaMediana")])


summary(data_summary[, c("Longitud", "Latitud")])


punto_referencia <- st_sfc(st_point(c(-75, 5)), crs = 4326)

candidatos <- c("MAGNA-SIRGAS / Colombia Bogota zone (3116)" = 3116,
                "MAGNA-SIRGAS 2018 / Origen-Nacional (9377)" = 9377,
                "WGS84 / UTM zona 18N (32618)"                = 32618)

purrr::map_dfr(candidatos, function(epsg) {
  co <- st_coordinates(st_transform(punto_referencia, epsg))
  tibble(este = co[1], norte = co[2])
}, .id = "sistema")



datos_sf <- st_as_sf(data_geo_tem, coords = c("Longitud", "Latitud"), crs = 4326, remove = FALSE)

datos_sf







library(geodata)
library(sf)
library(ggplot2)

# 1. Descargar departamentos de Colombia (nivel 1)
depto_colombia <- gadm(country = "COL", level = 1, path = tempdir()) %>% 
  st_as_sf() %>% 
  st_transform(crs = 4326)

# 2. Graficar mapa por departamentos con puntos en color continuo
ggplot() +
  geom_sf(data = depto_colombia, fill = "white", color = "gray70") +
  geom_sf(data = datos_sf, aes(color = TemperaturaMediana), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "inferno", name = "Temp. (°C)") +
  coord_sf(crs = 4326) +
  theme_bw() +
  labs(
    title = "Ubicación de Estaciones por Departamento",
    subtitle = "Color continuo: Temperatura (°C)",
    x = "Longitud",
    y = "Latitud"
  )


vc_temp <- variogram(TemperaturaMediana ~ 1, datos_sf, cloud = TRUE)

ggplot(as.data.frame(vc_temp), aes(x = dist, y = gamma)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred") +
  theme_test(base_size = 13) +
  labs(x = "Distancia geográfica (m)", y = expression((Z(s[i]) - Z(s[j]))^2),
       title = "Semivariograma cloud — temperatura") +
  theme(axis.title = element_text(face = "bold"))



geodist  <- dist(st_coordinates(datos_sf))
dist_temp <- dist(data_summary$TemperaturaMediana)

set.seed(12345)
mantel_vegan <- vegan::mantel(geodist, dist_temp, permutations = 999)
mantel_vegan




mantel_ade4 <- ade4::mantel.rtest(geodist, dist_temp, nrepet = 999)
mantel_ade4




distancia_maxima <- max(dist(st_coordinates(datos_sf))) / 2
ancho_intervalo  <- distancia_maxima / 12   # ~12 intervalos, balance resolución/estabilidad

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









