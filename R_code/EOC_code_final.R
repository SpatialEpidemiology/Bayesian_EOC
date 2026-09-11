#packages
library(tidyverse)
library(tidycensus)
library(sociome)
library(sf)
library(spdep)
library(cowplot)
library(gridGraphics)
library(tmap)
library(INLA)
library(inlatools)

tmap_mode("plot")


#select cancer and spatial files
#cancer_site must be "allsite", "crc", "breast", or "lung"
cancer_site <- "allsite"
early_file <- "EARLY CANCER DEATH SPATIAL FILE"
older_file <- "OLDER CANCER DEATH SPATIAL FILE"


#county shapefile
us <- unique(fips_codes$state)[1:51]
cou <- get_acs(year=2010, geography="county", state=us, survey="acs5",
               variables="B01001_001", geometry=TRUE)
cou$st <- substr(cou$GEOID, start=1, stop=2)
cou <- cou[which(cou$st != "78" & cou$st != "72" &
                   cou$st != "02" & cou$st != "15"), ]
cou <- st_transform(cou, crs=3857)
cou$one <- 1


#state shapefile
st <- get_acs(year=2014, geography="state", variables="B01001_001", geometry=TRUE)
st <- st_transform(st, crs=3857)
st$st <- substr(st$GEOID, start=1, stop=2)
st <- st[which(st$st != "72" & st$st != "78" &
                 st$st != "15" & st$st != "02"), ]
st$one <- 1


#read in data
rs_y_sf4 <- st_read(early_file, quiet=TRUE) %>%
  mutate(
    GEOID=str_pad(as.character(GEOID), 5, pad="0"),
    GEOID=if_else(GEOID == "46113", "46102", GEOID)
  )

rs_o_sf4 <- st_read(older_file, quiet=TRUE) %>%
  mutate(
    GEOID=str_pad(as.character(GEOID), 5, pad="0"),
    GEOID=if_else(GEOID == "46113", "46102", GEOID)
  )


#variables required for models
variables_y <- c("observed", "expected", "ADI", "b0049_pct", "h0049_pct",
                 "ob_cdc", "sm_mean", "al_mean", "unins_pct0054",
                 "LONGITUDE", "LATITUDE")
variables_o <- c("observed", "expected", "ADI", "b50_pct", "h50_pct",
                 "ob_cdc", "sm_mean", "al_mean", "unins_pct55",
                 "LONGITUDE", "LATITUDE")

if (cancer_site != "breast") {
  variables_y <- c(variables_y, "pct_fem")
  variables_o <- c(variables_o, "pct_fem")
}

if (!cancer_site %in% c("allsite", "crc", "breast", "lung")) {
  stop("cancer_site must be allsite, crc, breast, or lung")
}

missing_y <- setdiff(variables_y, names(rs_y_sf4))
missing_o <- setdiff(variables_o, names(rs_o_sf4))

if (length(missing_y) > 0) stop("Missing early cancer death variables: ", paste(missing_y, collapse=", "))
if (length(missing_o) > 0) stop("Missing older cancer death variables: ", paste(missing_o, collapse=", "))


#creating neighborhoods
#early cancer death
rs_y_sf <- st_drop_geometry(rs_y_sf4)
rs_y_sf <- st_as_sf(rs_y_sf, coords=c("LONGITUDE", "LATITUDE"), crs=4326)
knn <- knearneigh(st_coordinates(rs_y_sf), k=8)
nb <- knn2nb(knn)
nb <- make.sym.nb(nb)
nb2INLA(paste0("map_", cancer_site, "_early.adj"), nb)
gr <- inla.read.graph(paste0("map_", cancer_site, "_early.adj"))
rs_y_sf4$idarea <- seq_len(nrow(rs_y_sf4))

#older cancer death
rs_o_sf <- st_drop_geometry(rs_o_sf4)
rs_o_sf <- st_as_sf(rs_o_sf, coords=c("LONGITUDE", "LATITUDE"), crs=4326)
knn2 <- knearneigh(st_coordinates(rs_o_sf), k=8)
nb2 <- knn2nb(knn2)
nb2 <- make.sym.nb(nb2)
nb2INLA(paste0("map_", cancer_site, "_older.adj"), nb2)
gr2 <- inla.read.graph(paste0("map_", cancer_site, "_older.adj"))
rs_o_sf4$idarea <- seq_len(nrow(rs_o_sf4))


#priors set up based on rule of thumb (Simpson et al. 2017)
prior <- list(
  prec=list(
    prior="pc.prec",
    param=c(0.5 / 0.31, 0.01)),
  phi=list(
    prior="pc",
    param=c(0.5, 2 / 3))
)


#unadjusted formulas
formula_un_y <- observed ~ 1 +
  f(idarea, model="bym2", graph=gr, scale.model=TRUE, hyper=prior)

formula_un_o <- observed ~ 1 +
  f(idarea, model="bym2", graph=gr2, scale.model=TRUE, hyper=prior)


#adjusted formulas
#pct_fem is not included for breast cancer
if (cancer_site == "breast") {
  formula_y <- observed ~ ADI + b0049_pct + h0049_pct + ob_cdc + sm_mean +
    al_mean + unins_pct0054 +
    f(idarea, model="bym2", graph=gr, scale.model=TRUE, hyper=prior)

  formula_o <- observed ~ ADI + b50_pct + h50_pct + ob_cdc + sm_mean +
    al_mean + unins_pct55 +
    f(idarea, model="bym2", graph=gr2, scale.model=TRUE, hyper=prior)
} else {
  formula_y <- observed ~ ADI + b0049_pct + h0049_pct + ob_cdc + sm_mean +
    al_mean + unins_pct0054 + pct_fem +
    f(idarea, model="bym2", graph=gr, scale.model=TRUE, hyper=prior)

  formula_o <- observed ~ ADI + b50_pct + h50_pct + ob_cdc + sm_mean +
    al_mean + unins_pct55 + pct_fem +
    f(idarea, model="bym2", graph=gr2, scale.model=TRUE, hyper=prior)
}


#unadjusted models
#Poisson is used for the unadjusted models
#early cancer death
res_un <- inla(
  formula_un_y,
  family="poisson",
  data=st_drop_geometry(rs_y_sf4),
  E=rs_y_sf4$expected,
  control.predictor=list(compute=TRUE),
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE, cpo=TRUE),
  control.inla=list(int.strategy="ccd")
)

#older cancer death
res_un2 <- inla(
  formula_un_o,
  family="poisson",
  data=st_drop_geometry(rs_o_sf4),
  E=rs_o_sf4$expected,
  control.predictor=list(compute=TRUE),
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE, cpo=TRUE),
  control.inla=list(int.strategy="ccd")
)


#adjusted models
#negative binomial is used for the adjusted models
#early cancer death
res <- inla(
  formula_y,
  family="nbinomial",
  data=st_drop_geometry(rs_y_sf4),
  E=rs_y_sf4$expected,
  control.predictor=list(compute=TRUE),
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE,
                       cpo=TRUE, config=TRUE),
  control.inla=list(int.strategy="ccd")
)

#older cancer death
res2 <- inla(
  formula_o,
  family="nbinomial",
  data=st_drop_geometry(rs_o_sf4),
  E=rs_o_sf4$expected,
  control.predictor=list(compute=TRUE),
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE,
                       cpo=TRUE, config=TRUE),
  control.inla=list(int.strategy="ccd")
)


#model checks
res$dic$dic
res$waic$waic
res2$dic$dic
res2$waic$waic


#testing marginal plots
#early cancer death
marginal <- data.frame(inla.smarginal(res$marginals.fixed$unins_pct0054))
ggplot(marginal, aes(x=x, y=y)) +
  geom_line() +
  labs(x=expression(beta[1]), y="Density") +
  ggtitle("Posterior distribution of coefficient for uninsured") +
  geom_vline(xintercept=0, col="black") +
  theme_bw()

#older cancer death
marginal2 <- data.frame(inla.smarginal(res2$marginals.fixed$unins_pct55))
ggplot(marginal2, aes(x=x, y=y)) +
  geom_line() +
  labs(x=expression(beta[1]), y="Density") +
  ggtitle("Posterior distribution of coefficient for uninsured") +
  geom_vline(xintercept=0, col="black") +
  theme_bw()


#creating relative risks and exceedance probabilities
#early cancer death
rs_y_sf4$RR <- map_dbl(
  res$marginals.linear.predictor,
  ~inla.emarginal(exp, .x)
)
rs_y_sf4$LL <- map_dbl(
  res$marginals.linear.predictor,
  ~exp(inla.qmarginal(0.025, .x))
)
rs_y_sf4$UL <- map_dbl(
  res$marginals.linear.predictor,
  ~exp(inla.qmarginal(0.975, .x))
)
rs_y_sf4$exc <- map_dbl(
  res$marginals.linear.predictor,
  ~1 - inla.pmarginal(log(1.5), .x)
)
rs_y_sf4 <- st_transform(rs_y_sf4, crs=3857)

#older cancer death
rs_o_sf4$RR <- map_dbl(
  res2$marginals.linear.predictor,
  ~inla.emarginal(exp, .x)
)
rs_o_sf4$LL <- map_dbl(
  res2$marginals.linear.predictor,
  ~exp(inla.qmarginal(0.025, .x))
)
rs_o_sf4$UL <- map_dbl(
  res2$marginals.linear.predictor,
  ~exp(inla.qmarginal(0.975, .x))
)
rs_o_sf4$exc <- map_dbl(
  res2$marginals.linear.predictor,
  ~1 - inla.pmarginal(log(1.5), .x)
)
rs_o_sf4 <- st_transform(rs_o_sf4, crs=3857)


#RR plots
legend_title2 <- "RR"

#early cancer death
psup_RR_eo <- tm_shape(cou) +
  tm_polygons("one", legend.show=FALSE, alpha=0.5) +
  tm_shape(rs_y_sf4) +
  tm_borders(lwd=0) +
  tm_polygons("RR", palette="viridis",
              breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, Inf),
              legend.show=TRUE, title=legend_title2, colorNA="lightgrey") +
  tm_compass() +
  tm_scale_bar(position=0.78, width=0.1) +
  tm_credits("Pale yellow = Suppressed", position=c("center", "bottom")) +
  tm_layout(frame=FALSE)
psup_RR_eo

#older cancer death
psup_RR_ao <- tm_shape(cou) +
  tm_polygons("one", legend.show=FALSE, alpha=0.5) +
  tm_shape(rs_o_sf4) +
  tm_borders(lwd=0) +
  tm_polygons("RR", palette="viridis",
              breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, Inf),
              legend.show=TRUE, title=legend_title2, colorNA="lightgrey") +
  tm_compass() +
  tm_scale_bar(position=0.78, width=0.1) +
  tm_credits("Pale yellow = Suppressed", position=c("center", "bottom")) +
  tm_layout(frame=FALSE)
psup_RR_ao


#exceedance probability plots for RR>1.5
#early cancer death
plot(rs_y_sf4["exc"], border=NA, main=NA, key.pos=1,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
plot(cou["one"], border="darkgray", col="gray", main=NA,
     key.pos=NA, reset=FALSE, add=TRUE)
plot(rs_y_sf4["exc"], border=NA, main=NA, key.pos=NA,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
mtext("Early cancer death exceedance probability for RR > 1.5", side=1)
p_exc_e <- recordPlot()
p_exc_e

#older cancer death
plot(rs_o_sf4["exc"], border=NA, main=NA, key.pos=1,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
plot(cou["one"], border="darkgray", col="gray", main=NA,
     key.pos=NA, reset=FALSE, add=TRUE)
plot(rs_o_sf4["exc"], border=NA, main=NA, key.pos=NA,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
mtext("Older cancer death exceedance probability for RR > 1.5", side=1)
p_exc_o <- recordPlot()
p_exc_o


#RR 1.25 sensitivity analysis for early cancer death for all cancers
rs_y_sf4$exc125 <- map_dbl(
  res$marginals.linear.predictor,
  ~1 - inla.pmarginal(log(1.25), .x)
)

plot(rs_y_sf4["exc125"], border=NA, main=NA, key.pos=1,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
plot(cou["one"], border="darkgray", col="gray", main=NA,
     key.pos=NA, reset=FALSE, add=TRUE)
plot(rs_y_sf4["exc125"], border=NA, main=NA, key.pos=NA,
     pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
mtext("Early cancer death exceedance probability for RR > 1.25", side=1)
p_exc_e_125 <- recordPlot()
p_exc_e_125


#additional exceedance threshold sensitivity analyses are for all site cancer only
if (cancer_site == "allsite") {
  rs_y_sf4$exc100 <- map_dbl(
    res$marginals.linear.predictor,
    ~1 - inla.pmarginal(log(1), .x)
  )
  rs_y_sf4$exc175 <- map_dbl(
    res$marginals.linear.predictor,
    ~1 - inla.pmarginal(log(1.75), .x)
  )
  rs_y_sf4$exc200 <- map_dbl(
    res$marginals.linear.predictor,
    ~1 - inla.pmarginal(log(2), .x)
  )
  rs_y_sf4$exc225 <- map_dbl(
    res$marginals.linear.predictor,
    ~1 - inla.pmarginal(log(2.25), .x)
  )

  #all site early cancer death at RR>1
  plot(rs_y_sf4["exc100"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", main=NA,
       key.pos=NA, reset=FALSE, add=TRUE)
  plot(rs_y_sf4["exc100"], border=NA, main=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  mtext("All site early cancer death: RR > 1", side=1)
  p_threshold_100 <- recordPlot()

  #all site early cancer death at RR>1.75
  plot(rs_y_sf4["exc175"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", main=NA,
       key.pos=NA, reset=FALSE, add=TRUE)
  plot(rs_y_sf4["exc175"], border=NA, main=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  mtext("All site early cancer death: RR > 1.75", side=1)
  p_threshold_175 <- recordPlot()

  #all site early cancer death at RR>2
  plot(rs_y_sf4["exc200"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", main=NA,
       key.pos=NA, reset=FALSE, add=TRUE)
  plot(rs_y_sf4["exc200"], border=NA, main=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  mtext("All site early cancer death: RR > 2", side=1)
  p_threshold_200 <- recordPlot()

  #all site early cancer death at RR>2.25
  plot(rs_y_sf4["exc225"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", main=NA,
       key.pos=NA, reset=FALSE, add=TRUE)
  plot(rs_y_sf4["exc225"], border=NA, main=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  mtext("All site early cancer death: RR > 2.25", side=1)
  p_threshold_225 <- recordPlot()

  plot_grid(
    p_threshold_100,
    p_threshold_175,
    p_threshold_200,
    p_threshold_225,
    nrow=2,
    labels=c("A", "B", "C", "D")
  )
}


#covariate year sensitivity analyses for all site cancer only
if (cancer_site == "allsite") {
  #ACS insurance variables
  insured_variables <- c(
    m_u6="B27001_004", m617="B27001_007", m1824="B27001_010",
    m2534="B27001_013", m3544="B27001_016", m4554="B27001_019",
    m5564="B27001_022", m6574="B27001_025", m75="B27001_028",
    f_u6="B27001_032", f617="B27001_035", f1824="B27001_038",
    f2534="B27001_041", f3544="B27001_044", f4554="B27001_047",
    f5564="B27001_050", f6574="B27001_053", f75="B27001_056"
  )

  uninsured_variables <- c(
    m_u6="B27001_005", m617="B27001_008", m1824="B27001_011",
    m2534="B27001_014", m3544="B27001_017", m4554="B27001_020",
    m5564="B27001_023", m6574="B27001_026", m75="B27001_029",
    f_u6="B27001_033", f617="B27001_036", f1824="B27001_039",
    f2534="B27001_042", f3544="B27001_045", f4554="B27001_048",
    f5564="B27001_051", f6574="B27001_054", f75="B27001_057"
  )

  #2012-2016 ACS covariates
  adi_2016 <- get_adi(
    geography="county", state=us, dataset="acs5", year=2016
  ) %>%
    transmute(
      GEOID=str_pad(as.character(GEOID), 5, pad="0"),
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      ADI
    )

  demographic_2016 <- get_acs(
    year=2016, geography="county", state=us, survey="acs5",
    variables=c(
      race_total="B02001_001", black="B02001_003",
      ethnicity_total="B03003_001", hispanic="B03003_003",
      sex_total="B01001_001", female="B01001_026"
    )
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    transmute(
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      black_pct=100 * black / race_total,
      hispanic_pct=100 * hispanic / ethnicity_total,
      pct_fem=100 * female / sex_total
    )

  insured_2016 <- get_acs(
    year=2016, geography="county", state=us, survey="acs5",
    variables=insured_variables
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    mutate(insured=rowSums(pick(-GEOID))) %>%
    select(GEOID, insured)

  uninsured_2016 <- get_acs(
    year=2016, geography="county", state=us, survey="acs5",
    variables=uninsured_variables
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    mutate(uninsured=rowSums(pick(-GEOID))) %>%
    select(GEOID, uninsured)

  insurance_2016 <- insured_2016 %>%
    left_join(uninsured_2016, by="GEOID") %>%
    transmute(
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      unins_pct=100 * uninsured / (insured + uninsured)
    )

  covariates_2016 <- adi_2016 %>%
    full_join(demographic_2016, by="GEOID") %>%
    full_join(insurance_2016, by="GEOID")

  #2016-2020 ACS covariates
  adi_2020 <- get_adi(
    geography="county", state=us, dataset="acs5", year=2020
  ) %>%
    transmute(
      GEOID=str_pad(as.character(GEOID), 5, pad="0"),
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      ADI
    )

  demographic_2020 <- get_acs(
    year=2020, geography="county", state=us, survey="acs5",
    variables=c(
      race_total="B02001_001", black="B02001_003",
      ethnicity_total="B03003_001", hispanic="B03003_003",
      sex_total="B01001_001", female="B01001_026"
    )
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    transmute(
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      black_pct=100 * black / race_total,
      hispanic_pct=100 * hispanic / ethnicity_total,
      pct_fem=100 * female / sex_total
    )

  insured_2020 <- get_acs(
    year=2020, geography="county", state=us, survey="acs5",
    variables=insured_variables
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    mutate(insured=rowSums(pick(-GEOID))) %>%
    select(GEOID, insured)

  uninsured_2020 <- get_acs(
    year=2020, geography="county", state=us, survey="acs5",
    variables=uninsured_variables
  ) %>%
    select(GEOID, variable, estimate) %>%
    pivot_wider(names_from=variable, values_from=estimate) %>%
    mutate(uninsured=rowSums(pick(-GEOID))) %>%
    select(GEOID, uninsured)

  insurance_2020 <- insured_2020 %>%
    left_join(uninsured_2020, by="GEOID") %>%
    transmute(
      GEOID=if_else(GEOID == "46113", "46102", GEOID),
      unins_pct=100 * uninsured / (insured + uninsured)
    )

  covariates_2020 <- adi_2020 %>%
    full_join(demographic_2020, by="GEOID") %>%
    full_join(insurance_2020, by="GEOID")

  #replace primary ACS covariates, retaining outcomes and fixed covariates
  rs_y_2016 <- rs_y_sf4 %>%
    select(-any_of(c("ADI", "b0049_pct", "h0049_pct",
                     "unins_pct0054", "pct_fem", "idarea"))) %>%
    left_join(covariates_2016, by="GEOID")

  rs_o_2016 <- rs_o_sf4 %>%
    select(-any_of(c("ADI", "b50_pct", "h50_pct",
                     "unins_pct55", "pct_fem", "idarea"))) %>%
    left_join(covariates_2016, by="GEOID")

  rs_y_2020 <- rs_y_sf4 %>%
    select(-any_of(c("ADI", "b0049_pct", "h0049_pct",
                     "unins_pct0054", "pct_fem", "idarea"))) %>%
    left_join(covariates_2020, by="GEOID")

  rs_o_2020 <- rs_o_sf4 %>%
    select(-any_of(c("ADI", "b50_pct", "h50_pct",
                     "unins_pct55", "pct_fem", "idarea"))) %>%
    left_join(covariates_2020, by="GEOID")

  #2012-2016 for early cancer death
  points_y_2016 <- st_drop_geometry(rs_y_2016)
  points_y_2016 <- st_as_sf(points_y_2016,
                            coords=c("LONGITUDE", "LATITUDE"), crs=4326)
  knn_y_2016 <- knearneigh(st_coordinates(points_y_2016), k=8)
  nb_y_2016 <- make.sym.nb(knn2nb(knn_y_2016))
  nb2INLA("map_allsite_2016_early.adj", nb_y_2016)
  gr_y_2016 <- inla.read.graph("map_allsite_2016_early.adj")
  rs_y_2016$idarea <- seq_len(nrow(rs_y_2016))

  #2012-2016 for older cancer death
  points_o_2016 <- st_drop_geometry(rs_o_2016)
  points_o_2016 <- st_as_sf(points_o_2016,
                            coords=c("LONGITUDE", "LATITUDE"), crs=4326)
  knn_o_2016 <- knearneigh(st_coordinates(points_o_2016), k=8)
  nb_o_2016 <- make.sym.nb(knn2nb(knn_o_2016))
  nb2INLA("map_allsite_2016_older.adj", nb_o_2016)
  gr_o_2016 <- inla.read.graph("map_allsite_2016_older.adj")
  rs_o_2016$idarea <- seq_len(nrow(rs_o_2016))

  #2016-2020 for early cancer death
  points_y_2020 <- st_drop_geometry(rs_y_2020)
  points_y_2020 <- st_as_sf(points_y_2020,
                            coords=c("LONGITUDE", "LATITUDE"), crs=4326)
  knn_y_2020 <- knearneigh(st_coordinates(points_y_2020), k=8)
  nb_y_2020 <- make.sym.nb(knn2nb(knn_y_2020))
  nb2INLA("map_allsite_2020_early.adj", nb_y_2020)
  gr_y_2020 <- inla.read.graph("map_allsite_2020_early.adj")
  rs_y_2020$idarea <- seq_len(nrow(rs_y_2020))

  #2016-2020 for older cancer death
  points_o_2020 <- st_drop_geometry(rs_o_2020)
  points_o_2020 <- st_as_sf(points_o_2020,
                            coords=c("LONGITUDE", "LATITUDE"), crs=4326)
  knn_o_2020 <- knearneigh(st_coordinates(points_o_2020), k=8)
  nb_o_2020 <- make.sym.nb(knn2nb(knn_o_2020))
  nb2INLA("map_allsite_2020_older.adj", nb_o_2020)
  gr_o_2020 <- inla.read.graph("map_allsite_2020_older.adj")
  rs_o_2020$idarea <- seq_len(nrow(rs_o_2020))

  #sensitivity model formulas
  formula_y_2016 <- observed ~ ADI + black_pct + hispanic_pct + ob_cdc +
    sm_mean + al_mean + unins_pct + pct_fem +
    f(idarea, model="bym2", graph=gr_y_2016,
      scale.model=TRUE, hyper=prior)
  formula_o_2016 <- observed ~ ADI + black_pct + hispanic_pct + ob_cdc +
    sm_mean + al_mean + unins_pct + pct_fem +
    f(idarea, model="bym2", graph=gr_o_2016,
      scale.model=TRUE, hyper=prior)
  formula_y_2020 <- observed ~ ADI + black_pct + hispanic_pct + ob_cdc +
    sm_mean + al_mean + unins_pct + pct_fem +
    f(idarea, model="bym2", graph=gr_y_2020,
      scale.model=TRUE, hyper=prior)
  formula_o_2020 <- observed ~ ADI + black_pct + hispanic_pct + ob_cdc +
    sm_mean + al_mean + unins_pct + pct_fem +
    f(idarea, model="bym2", graph=gr_o_2020,
      scale.model=TRUE, hyper=prior)

  #sensitivity models
  res_y_2016 <- inla(
    formula_y_2016, family="nbinomial",
    data=st_drop_geometry(rs_y_2016), E=rs_y_2016$expected,
    control.predictor=list(compute=TRUE),
    control.compute=list(return.marginals.predictor=TRUE,
                         dic=TRUE, waic=TRUE, cpo=TRUE),
    control.inla=list(int.strategy="ccd")
  )
  res_o_2016 <- inla(
    formula_o_2016, family="nbinomial",
    data=st_drop_geometry(rs_o_2016), E=rs_o_2016$expected,
    control.predictor=list(compute=TRUE),
    control.compute=list(return.marginals.predictor=TRUE,
                         dic=TRUE, waic=TRUE, cpo=TRUE),
    control.inla=list(int.strategy="ccd")
  )
  res_y_2020 <- inla(
    formula_y_2020, family="nbinomial",
    data=st_drop_geometry(rs_y_2020), E=rs_y_2020$expected,
    control.predictor=list(compute=TRUE),
    control.compute=list(return.marginals.predictor=TRUE,
                         dic=TRUE, waic=TRUE, cpo=TRUE),
    control.inla=list(int.strategy="ccd")
  )
  res_o_2020 <- inla(
    formula_o_2020, family="nbinomial",
    data=st_drop_geometry(rs_o_2020), E=rs_o_2020$expected,
    control.predictor=list(compute=TRUE),
    control.compute=list(return.marginals.predictor=TRUE,
                         dic=TRUE, waic=TRUE, cpo=TRUE),
    control.inla=list(int.strategy="ccd")
  )

  #sensitivity exceedance probabilities
  rs_y_2016$exc125 <- map_dbl(res_y_2016$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.25), .x))
  rs_y_2016$exc150 <- map_dbl(res_y_2016$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.5), .x))
  rs_o_2016$exc125 <- map_dbl(res_o_2016$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.25), .x))
  rs_o_2016$exc150 <- map_dbl(res_o_2016$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.5), .x))
  rs_y_2020$exc125 <- map_dbl(res_y_2020$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.25), .x))
  rs_y_2020$exc150 <- map_dbl(res_y_2020$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.5), .x))
  rs_o_2020$exc125 <- map_dbl(res_o_2020$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.25), .x))
  rs_o_2020$exc150 <- map_dbl(res_o_2020$marginals.linear.predictor,
                              ~1 - inla.pmarginal(log(1.5), .x))

  rs_y_2016 <- st_transform(rs_y_2016, crs=3857)
  rs_o_2016 <- st_transform(rs_o_2016, crs=3857)
  rs_y_2020 <- st_transform(rs_y_2020, crs=3857)
  rs_o_2020 <- st_transform(rs_o_2020, crs=3857)

  #2012-2016 early cancer death at RR>1.25
  plot(rs_y_2016["exc125"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_y_2016["exc125"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_y_2016_125 <- recordPlot()

  #2012-2016 early cancer death at RR>1.5
  plot(rs_y_2016["exc150"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_y_2016["exc150"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_y_2016_150 <- recordPlot()

  #2016-2020 early cancer death at RR>1.25
  plot(rs_y_2020["exc125"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_y_2020["exc125"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_y_2020_125 <- recordPlot()

  #2016-2020 early cancer death at RR>1.5
  plot(rs_y_2020["exc150"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_y_2020["exc150"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_y_2020_150 <- recordPlot()

  #all site covariate year sensitivity figure for early cancer death
  plot_grid(p_y_2016_125, p_y_2016_150, p_y_2020_125, p_y_2020_150,
            nrow=2, labels=c("A", "B", "C", "D"))

  #2012-2016 older cancer death at RR>1.25
  plot(rs_o_2016["exc125"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_o_2016["exc125"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_o_2016_125 <- recordPlot()

  #2012-2016 older cancer death at RR>1.5
  plot(rs_o_2016["exc150"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_o_2016["exc150"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_o_2016_150 <- recordPlot()

  #2016-2020 older cancer death at RR>1.25
  plot(rs_o_2020["exc125"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_o_2020["exc125"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_o_2020_125 <- recordPlot()

  #2016-2020 older cancer death at RR>1.5
  plot(rs_o_2020["exc150"], border=NA, main=NA, key.pos=1,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE)
  plot(cou["one"], border="darkgray", col="gray", key.pos=NA,
       reset=FALSE, add=TRUE)
  plot(rs_o_2020["exc150"], border=NA, key.pos=NA,
       pal=hcl.colors(10, "Heat", rev=TRUE), reset=FALSE, add=TRUE)
  plot(st["one"], border="white", col=sf.colors(1, alpha=0), add=TRUE)
  p_o_2020_150 <- recordPlot()

  #all site covariate year sensitivity figure for older cancer death
  plot_grid(p_o_2016_125, p_o_2016_150, p_o_2020_125, p_o_2020_150,
            nrow=2, labels=c("A", "B", "C", "D"))
}
