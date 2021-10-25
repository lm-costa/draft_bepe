Draft
================

# **Introduction**

The data are from the Fluox, and are already processed for the clear sky
days, also the database are filtered for the DOY&gt;=128, the bellow
code was apply to be abble of uploading this data to github because the
size of the first file was larger than 100Mb:

``` r
# 
#namefile <- list.files("data-raw/spectra/", pattern="JB-009-ESA_F")
# 
# ori_tab <- readr::read_csv(paste0("data-raw/spectra/",namefile),col_names = TRUE)
# col_Names <- as.character(ori_tab$F_final1)
# 
# trans_tab <- t(ori_tab)
# colnames(trans_tab) <- col_Names
# trans_tab <- as.data.frame(trans_tab)
# 
# trans_tab <- trans_tab[-1, ]
# trans_tab <- trans_tab[,-2]
# colnames(trans_tab)[1] <- "DOY"
# 
# trans_tab <- trans_tab |>
#   dplyr::filter(
#     DOY > 127.9
#   )
# namefile_rds <- sub(".csv",".rds",namefile)
# readr::write_rds(trans_tab, paste0("data/",namefile_rds))
```

Now we are keeping going from here with the `.rds` file

``` r
namefile <- list.files("data-raw/", pattern = "_F_")
ori_tab <- readr::read_rds(paste0("data-raw/",namefile))
```

Reading the Index file

``` r
ind_file <- readr::read_csv("data-raw/table_all_index.csv") |>
  dplyr::filter(DOYdayfrac > 127.9)
```

# **Processing**

First we are gonna make the time table and made a mean for every 30 m

``` r
time_table <- cbind(ind_file,ori_tab)

time_table <- time_table |>
  dplyr::filter(lubridate::hour(time_table$UTC_datetime) > 8)

min_int<-"30 min" #step 
sttime <- time_table$UTC_datetime[1] # start time of the serie
endtime <- time_table$UTC_datetime[13577] #final time of the serie


## Making the agreggation

timetoagg<-seq.POSIXt(from = as.POSIXct(sttime,tz = "UTC"),
                      to = as.POSIXct(endtime,tz = "UTC"),
                      by = min_int)

time_agg <- aggregate(x = time_table,
                      by = list(cut(as.POSIXct(time_table$UTC_datetime,tz = "UTC"),
                                    timetoagg)), FUN = mean)
```

Now we need to load the GPP data and find the days and hours that bolth
have measurments

``` r
gppFile <- list.files("data-raw/", pattern = "aa_")
gpp_table <- readr::read_csv(paste0("data-raw/", gppFile))
```

If we try to use the function `dplyr::semi_join()` we are not gonna get
any result, because the data frames don’t have any column to compare, so
we must create in this case

``` r
time_agg <- time_agg |> dplyr::mutate(day = time_agg$DOY%/%1,
                                      Hr1 = lubridate::hour(as.character(Group.1)),
                                      Mn1 = lubridate::minute(as.character(Group.1)))

time_agg <- time_agg[,-7] # this line is because they have the same same in the GPP table, but with no match 

gpp_table <- gpp_table |> dplyr::mutate(day=DOY%/%1)

# Now we can use the semi_join

comp_f_gpp <- dplyr::semi_join(time_agg,gpp_table)
comp_gpp_f <- dplyr::semi_join(gpp_table,time_agg)
comp_gpp_f <- comp_gpp_f[,-c(6,176)]

df_final <- cbind(comp_f_gpp,comp_gpp_f)
```

Now we are gonna do a integral in the F spectra, but first we have to do
some transformations in the `time_agg`.

``` r
spectra_Table <- time_agg[,-c(689:691)]
spectra_Table <- as.data.frame(t(spectra_Table))
c_spectra <- paste0(lubridate::hour(as.character(spectra_Table[1,])),":",
                                    lubridate::minute(as.character(spectra_Table[1, ])))
colnames(spectra_Table) <- c_spectra
spectra_Table <- spectra_Table[-c(1:6), ]
rownames(spectra_Table) <- seq(1:682)
spectra <- names(ori_tab[0, c(2:683)])
spectra_Table <- cbind(spectra,spectra_Table)
colnames(spectra_Table) <- make.unique(names(spectra_Table))
spectra_Table[,1] <-as.numeric(spectra_Table[,1])
```

Integrating function:

``` r
integral <- function(x,y){
  approxfun(x,y)
}

result <- vector("list", length = ncol(spectra_Table))

for (i in 1:(ncol(spectra_Table))){
  int <- integral(spectra_Table$spectra, spectra_Table[,i])
  result[[i]] <- integrate(int, min(spectra_Table$spectra),
                           max(spectra_Table$spectra))
}

for(i in 1:length(result)){
  if(i==1){
    fint <- result[[i]]$value
  }else{
    fint_a <- result[[i]]$value
    fint <- rbind(fint,fint_a)
    fint <- as.data.frame(fint)
  }
}

row.names(fint) <- seq(1:377)
fint <- as.data.frame(fint[-1,]) # the first observasation is the integral of wavelenght in function of it self

colnames(fint) <- "Fint"

## Now we are gonna add the integral in the final dataset

df_final <- cbind(df_final,fint)
```

Before going plot the data, first we are gonna make a subset for some
variables

``` r
df_final_sub <- df_final |>
  dplyr::select(
    Group.1, day, Hr1, L750_Wm2nm1sr1, SZA, PAR_Wm2,
    `685.427527917867`,`743.81631095515`, Fint,
    FCO2,ET,NEE,Tair,Tsoil,NEE_orig,rH_orig,
    GPP_f, GPP_U05_f, GPP_U50_f,GPP_U95_f, GPP_DT, GPP_DT_SD,
    GPP_DT_U95, GPP_DT_U95_SD, GPP_DT_U50, GPP_DT_U50_SD,GPP_DT_U05, GPP_DT_U05_SD

  ) |> 
  dplyr::glimpse()
```

    ## Rows: 376
    ## Columns: 28
    ## $ Group.1            <fct> 2018-05-12 09:00:12, 2018-05-12 09:30:12, 2018-05-1~
    ## $ day                <dbl> 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 1~
    ## $ Hr1                <int> 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 1~
    ## $ L750_Wm2nm1sr1     <dbl> 0.30742885, 0.32451641, 0.33455608, 0.34120232, 0.3~
    ## $ SZA                <dbl> 34.95596, 30.81628, 27.45133, 25.35580, 24.71953, 2~
    ## $ PAR_Wm2            <dbl> 364.37724, 387.53816, 403.78000, 413.90729, 417.906~
    ## $ `685.427527917867` <dbl> 1.0169006, 1.0777418, 1.0466246, 1.0804593, 1.08365~
    ## $ `743.81631095515`  <dbl> 4.089323, 4.242339, 4.425685, 4.502813, 4.481843, 4~
    ## $ Fint               <dbl> 228.78444, 238.46321, 245.18873, 250.34433, 249.735~
    ## $ FCO2               <dbl> -7.08e-07, -6.87e-07, -6.93e-07, -5.58e-07, -5.08e-~
    ## $ ET                 <dbl> NA, NA, 0.000122369, NA, NA, NA, 0.000097600, NA, N~
    ## $ NEE                <dbl> -16.1000, -15.6000, -15.7000, -12.7000, -11.6000, -~
    ## $ Tair               <dbl> NA, NA, 19.9099, NA, NA, NA, 22.9810, NA, NA, NA, N~
    ## $ Tsoil              <dbl> NA, NA, 21.8680, NA, NA, NA, 24.3686, NA, NA, NA, N~
    ## $ NEE_orig           <dbl> -16.0820, -15.6160, -15.7490, -12.6750, -11.5520, -~
    ## $ rH_orig            <dbl> 82.8187, 81.1107, 73.8973, 72.4677, 70.3815, 66.623~
    ## $ GPP_f              <dbl> 24.7, 24.6, 24.5, 22.0, 21.3, 23.2, 20.3, 18.2, 20.~
    ## $ GPP_U05_f          <dbl> 23.7, 23.6, 23.6, 21.0, 20.2, 22.1, 19.2, 17.2, 19.~
    ## $ GPP_U50_f          <dbl> 24.6, 24.4, 24.4, 21.8, 21.1, 23.0, 20.1, 18.0, 20.~
    ## $ GPP_U95_f          <dbl> 26.1, 26.1, 26.0, 23.7, 23.1, 25.1, 22.1, 19.6, 22.~
    ## $ GPP_DT             <dbl> 18.8, 19.5, 20.1, 20.5, 20.5, 20.5, 20.8, 19.3, 20.~
    ## $ GPP_DT_SD          <dbl> 1.4265, 1.4120, 1.3955, 1.3790, 1.3648, 1.3516, 1.3~
    ## $ GPP_DT_U95         <dbl> 18.7, 19.4, 20.1, 20.5, 20.4, 20.5, 20.8, 19.2, 20.~
    ## $ GPP_DT_U95_SD      <dbl> 1.8851, 1.8744, 1.8599, 1.8431, 1.8264, 1.8102, 1.7~
    ## $ GPP_DT_U50         <dbl> 17.8, 18.5, 19.1, 19.5, 19.5, 19.6, 19.8, 18.3, 19.~
    ## $ GPP_DT_U50_SD      <dbl> 1.8059, 1.7984, 1.7868, 1.7726, 1.7575, 1.7427, 1.7~
    ## $ GPP_DT_U05         <dbl> 19.1, 19.8, 20.4, 20.7, 20.7, 20.8, 21.0, 19.6, 20.~
    ## $ GPP_DT_U05_SD      <dbl> 1.2953, 1.2757, 1.2553, 1.2366, 1.2223, 1.2097, 1.2~

Also we need another transformations to a correlation plot for all the
Fluoresce spectra

``` r
fint_t <- vector("numeric", length = 377)
fint_t[2:377] <- t(fint)

spectra_Table[683, ] <- fint_t
gpp <- vector("numeric", length = 377)
gpp[2:377] <- comp_gpp_f$GPP_DT_U95
spectra_Table[684, ] <- gpp

for(i in 2:ncol(spectra_Table)){
  if(i==2){
    tab <- as.numeric(spectra_Table[,i])
  }else{
    tab_a <- as.numeric(spectra_Table[,i])
    tab <- cbind(tab,tab_a)
    }
}
tab <- as.data.frame(tab)
colnames(tab) <- comp_gpp_f$TIMESTAMP
tab <- as.data.frame(t(tab))
colnames(tab) <- spectra_Table$spectra
colnames(tab)[683] <- "fint"
colnames(tab)[684] <- "gpp"

## 

pear_resul <- vector("list", length = ncol(tab))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab$gpp, tab[,i])
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl <- pear_resul[[i]]$estimate
  }else{
    correl_a <- pear_resul[[i]]$estimate
    correl <- rbind(correl,correl_a)
    correl <- as.data.frame(correl)
  }
}

rownames(correl) <- colnames(tab)
correl[,2] <- spectra_Table$spectra
```

# **Plots**

## **Correlation heat-map**

``` r
df_final_sub[,-1] |>
  cor(use = "p") |>
  corrplot::corrplot()
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Ploting correlation in the Fluoresce spectra, without F integral.
`F integral r = 0.65`

``` r
correl[-c(683:684),] |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(V2), y = cor))+
  ggplot2::geom_point()+
  ggplot2::xlab(label="Fluorescence Spectra in nm")
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## **Loop for se changes in the wave curve**

I just gonna do 4 plots, but you can apply this for any number

``` r
# starts from 2 because the first column is the spectra 

for(i in 2:4){
  print(spectra_Table|>
          dplyr::filter(spectra != 0) |> 
          ggplot2::ggplot(ggplot2::aes(x= spectra,
                                       y =as.numeric(spectra_Table[-c(683:684),i])
                                       ))+
          ggplot2::geom_point()+
          ggplot2::theme_bw())
}
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

## **Loop for the change in time for the wavelenght**

``` r
for(i in 7:9){
  print(time_agg |>
          ggplot2::ggplot(ggplot2::aes(x= UTC_datetime,
                                       y = time_agg[,i]))+
          ggplot2::geom_point()+
          ggplot2::geom_smooth()+
          ggplot2::theme_bw())
}
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

## **Regressions**

``` r
df_final_sub |> 
  ggplot2::ggplot(ggplot2::aes(x=GPP_DT_U95, y= Fint))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
  label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
df_final_sub |>
  dplyr::filter(Hr1 > 10 & Hr1 < 14) |> 
  ggplot2::ggplot(ggplot2::aes(x=GPP_DT_U95, y= Fint))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
  label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
df_final_sub |> 
  ggplot2::ggplot(ggplot2::aes(x=GPP_DT_U95, y= `743.81631095515`))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
  label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
df_final_sub |> 
  ggplot2::ggplot(ggplot2::aes(x=GPP_DT_U95, y= `685.427527917867`))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
  label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
df_final_sub |> 
  ggplot2::ggplot(ggplot2::aes(x=GPP_DT_U95, y= Fint))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggplot2::facet_wrap(~day)+
  ggplot2::theme_classic()+
  ggpubr::stat_regline_equation(ggplot2::aes(
  label =  paste(..rr.label..)))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

# **Clounds and reflectance**

Detection of clounds using the incoming radiance at 750

``` r
for(i in c(132,139,145,146,161,171,177,181,182,183,187,188,189,190,191,
           199,200,201,217,218,221,224,239,240)){
  print(time_agg |>
          dplyr::filter(DOYdayfrac%/%1 == i) |>
          ggplot2::ggplot(ggplot2::aes(x=DOYdayfrac, y= L750_Wm2nm1sr1))+
          ggplot2::geom_jitter()+
          ggplot2::geom_smooth())
  }
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-5.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-6.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-7.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-8.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-9.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-10.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-11.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-12.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-13.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-14.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-15.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-16.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-17.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-18.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-19.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-20.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-21.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-22.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-23.png)<!-- -->

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-20-24.png)<!-- -->

Here we are gonna see how the reflectance of the canopy are changing,
and for this we are gonna use the NDVI and NDVI red-edge

``` r
rfile <- list.files("data-raw/", pattern = "_R_")
rtab <- readr::read_rds(paste0("data-raw/",rfile))

rspec <- as.data.frame(t(rtab))
c_name <- rspec[1,]
colnames(rspec) <- c_name
rspec <- rspec[-1, ]
```

Using only the bands that are necessary

``` r
rnir <- rspec[c(616:682),]
rnir_spec <- as.numeric(rownames(rnir))
rnir <- cbind(rnir_spec,rnir)

rred <- rspec[c(1:49),]
rred_spec <- as.numeric(rownames(rred))
rred <- cbind(rred_spec,rred)

redge <- rspec[c(207:267),]
redge_spec <- as.numeric(rownames(redge))
redge <- cbind(redge_spec,redge)
```

doing the integral for these bands

``` r
rnir_result <- vector("list", length = ncol(rnir))
rred_result <- vector("list", length = ncol(rred))
redge_result <- vector("list", length = ncol(redge))

for (i in 1:(ncol(rnir))){
  int <- integral(rnir$rnir_spec, rnir[,i])
  rnir_result[[i]] <- integrate(int, min(rnir$rnir_spec), max(rnir$rnir_spec))
}

for(i in 1:length(rnir_result)){
  if(i==1){
    rnir_int <- rnir_result[[i]]$value
  }else{
    rnir_int_a <- rnir_result[[i]]$value
    rnir_int <- rbind(rnir_int,rnir_int_a)
    rnir_int <- as.data.frame(rnir_int)
  }
}
rnir_int <- rnir_int[-1,]
rnir_int <- as.data.frame(rnir_int)

for (i in 1:(ncol(rred))){
  int <- integral(rred$rred_spec, rred[,i])
  rred_result[[i]] <- integrate(int, min(rred$rred_spec), max(rred$rred_spec))
}

for(i in 1:length(rred_result)){
  if(i==1){
    rred_int <- rred_result[[i]]$value
  }else{
    rred_int_a <- rred_result[[i]]$value
    rred_int <- rbind(rred_int,rred_int_a)
    rred_int <- as.data.frame(rred_int)
  }
}
rred_int <- rred_int[-1,]
rred_int <- as.data.frame(rred_int)

for (i in 1:(ncol(redge))){
  int <- integral(redge$redge_spec, redge[,i])
  redge_result[[i]] <- integrate(int, min(redge$redge_spec), max(redge$redge_spec))
}

for(i in 1:length(redge_result)){
  if(i==1){
    redge_int <- redge_result[[i]]$value
  }else{
    redge_int_a <- redge_result[[i]]$value
    redge_int <- rbind(redge_int,redge_int_a)
    redge_int <- as.data.frame(redge_int)
  }
}

redge_int <- redge_int[-1,]
redge_int <- as.data.frame(redge_int)
```

Now the formula for the both ndvi

``` r
ndvi <- (rnir_int - rred_int)/(rnir_int+rred_int)
colnames(ndvi) <- "NDVI"

ndvi_r <- (rnir_int- redge_int)/(rnir_int + redge_int)
colnames(ndvi_r) <- "NDVI Red"

rtab <- cbind(rtab,ndvi,ndvi_r)

ndvi <- rtab |>
  dplyr::mutate(day = DOY%/%1) |>
  dplyr::group_by(day) |>
  dplyr::summarise(ndvi = mean(NDVI))

ndvi_r <- rtab |>
  dplyr::mutate(day = DOY%/%1) |>
  dplyr::group_by(day) |>
  dplyr::summarise(ndvi_r = mean(`NDVI Red`))
```

And Finaly the plots to see the changes

## **NDVI**

``` r
ndvi |>
  ggplot2::ggplot(ggplot2::aes(x=day, y=ndvi))+
  ggplot2::geom_jitter()+
  ggplot2::geom_point(data = subset(ndvi, ndvi< 0.9),
                      ggplot2::aes(x=day,y=ndvi)
                      , colour="red")
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## **NDVI Red Edge**

``` r
ndvi_r |>
  ggplot2::ggplot(ggplot2::aes(x=day, y=ndvi_r))+
  ggplot2::geom_jitter()+
  ggplot2::geom_point(data = subset(ndvi_r, ndvi_r< 0.59),
                      ggplot2::aes(x=day,y=ndvi_r)
                      , colour="red")
```

![](README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

# **Interpreting**

the following figure show to us, the days in the data-set that have
problems and the most provable cause, this days we are gonna exclud for
our analyses.

<img src="data/fig/Imagem1.png" width="1000px" style="display: block; margin: auto;" />
