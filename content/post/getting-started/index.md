---
title: Reconstruction life table functions from life expectancy

# Summary for listings and search engines
summary: This R code might help you reconstructing life table functions when only ex is known

# Link this post with a project
projects: []

# Date published
date: "2021-07-04T00:00:00Z"

# Date updated
lastmod: "2021-07-04T00:00:00Z"

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

# Featured image
# Place an image named `featured.jpg/png` in this page's folder and customize its options here.
image:
  caption: 'Plot'
  focal_point: ""
  placement: 2
  preview_only: false

authors:
- admin

tags:
- R coding

categories:
- Demography
---

## Deriving education-specific life tables with an iterative process
The life table survivors at age $x$ ($l_x$) can be obtained from life expectancy estimates at age $x$ ($e_x$) after assuming that in each age interval $x$ to $x+1$, people dying within this period live on average $1/2$ person-years ($a_x=0.5$):
\begin{equation}
l_{x+1}=\frac{l_x \cdot (2 \cdot e_x -1)}{1+2 \cdot e_{x+1}}.
\end{equation}
Please note, $l_0$ denotes the life table radix (usually defined as 100 000) and does not require estimation. Thus, the life table reconstruction starts with deriving $l_1$:
\begin{equation}
l_{1}=\frac{l_0 \cdot (2 \cdot e_0 -1)}{1+2 \cdot e_{1}}.
\end{equation}
In this way, the life table survivors at age 1 can be estimated from three known life table functions, i.e.,  $l_0$, $e_0$, and $e_1$. In the next step, $l_2$ is estimated from $l_1$, $e_1$, and $e_2$ and so forth. Once all $l_x$ are estimated on the basis of this algorithm, the remaining life table functions can be easily derived, such as $L_x$ ($L_x=(l_x+l_{x+1})/2)$. Theoretically, equation 1 enables us to reconstruct life table functions based on $e_x$ values (under the $a_x$ = 0.5 assumption). In practice, however, the reconstruction might require additional steps. For example, the $e_x$ values provided by Eurostat have only one decimal place. This limits the accuracy of the $l_x$ derivation and might result in constant $l_x$ values for several ages. To overcome this issue, we fitted a non-parametric curve to the data and predicted $e_x$ values with more decimal places. More specifically, we used the loess() function in R in order to obtain $e_x$ values with more decimal places that are as close as possible to the original $e_x$ values. In some cases, e.g., for the highly educated subpopulation in very low-mortality countries, the proposed derivation procedure still produces constant $l_x$ values at young ages. We solved this issue by focusing on $e_{30}$ and HLY at age 30.


The following code provides an example for calculating education-specific life tables when only the education-specific $e_x$ values are known. In other words, the aim of the code is to calculate the life table backwards, namely from $e_x$ to $p_x$. This is necessary because Eurostat does not provide education-specific life tables, but education-specific $e_x$ values are available. Please note, the results in this example will differ from the results in my paper (Sauerberg 2021) due to updates in the Eurostat database.
```r
library(dplyr)
library(eurostat)
#please load these packages and download the data like this:
data <- get_eurostat("demo_mlexpecedu", time_format = "num")

#rename and redefine the file 
data$isced11 <- as.character(data$isced11)
data$isced11 <- ifelse(data$isced11=="ED0-2", "lower", data$isced11)
data$isced11 <- ifelse(data$isced11=="ED3_4", "middle", data$isced11)
data$isced11 <- ifelse(data$isced11=="ED5-8", "higher", data$isced11)
data$isced11 <- ifelse(data$isced11=="TOTAL", "total", data$isced11)

data$age <- as.character(data$age)
data$age <- ifelse(data$age=="Y_LT1", "Y0", data$age)
data$age <- ifelse(data$age=="Y_GE85", "Y85", data$age)
data$age <- substring(data$age, 2)

data <- data[,-1]
colnames(data) <- c("sex","age","edu","country","year","ex")
data$age <- as.numeric(data$age)
#Filter for the year 2016, as we have done
data <- filter(data, year==2016)
```
The following function has the arguments "country.select", "edu.select" and "sex.select". Thus, the funcation allows to derive life tables for each educational level (high, middle, low, and total), for each country with available data (16 European countries), separated for men and women.

```r
my.function <- function(country.select, edu.select, sex.select) {

    select.country <- arrange(filter(data, country==country.select ,edu==edu.select &
                                               sex==sex.select),age)

#smooth to get more decimals by applying the loess function,
#then predict ex with more decimals
    grab.LE <- select.country$ex
    smooth.it <- loess(grab.LE~select.country$age, span=0.2)
    predict.it <- predict(smooth.it, seq(0,85,1))
    select.country$ex.decimals <- predict.it

    LT.derive <- data.frame(Age=0:85)
    LT.derive$lx <- NA

    LT.derive$ex <- select.country$ex.decimals
    LT.derive$lx[1] <- 100000
    LT.derive$Tx[1] <- 100000*select.country$ex.decimals[1]
    
#this loop refers to equation 1 in the paper
    for (j in 2:86) {

        upper <- LT.derive$lx[j-1]*(2*LT.derive$ex[j-1]-1)
        bottom <- 1+2*LT.derive$ex[j]
        LT.derive$lx[j] <- upper/bottom
    }
#Checks if lx is monotonic decreasing, i.e., no resurrection in the life table
    lx.diff <- diff(LT.derive$lx)
    lx.diff <- round(lx.diff, 5)

    if (all(diff(lx.diff) < 0)) {

        px <- c(LT.frame$lx[-1]/LT.frame$lx[-86],0)

    }else{
#sometimes, it is not, so I force it =)
#please note, this occurs usually at very young ages and won't affect
#LE at age 30 or older
        lx.diff[lx.diff>=0] <- -runif(length(lx.diff[lx.diff>=0]), 1, 5)
        lx.monotonic <- cumsum(c(100000, lx.diff))
        px <- c(lx.monotonic[-1]/lx.monotonic[-86],0)

        }
#from here, the life table is constructed very standardly
    lx <- round(c(100000, (cumprod(px)*100000)[1:(length(px)-1)]))
    dx <- round(c(-diff(lx), lx[length(lx)]))
    LT.derive$lx <- lx
    LT.derive$dx <- dx
    LT.derive$px <- px
    Lx1 <- lx[-1]+0.5[-length(px)]*dx[-length(dx)]
    Lx.open <- LT.derive$Tx[1]-sum(Lx1)
    LT.derive$Lx <- round(c(Lx1, Lx.open))
    LT.derive$Tx <- rev(cumsum(rev(LT.derive$Lx)))
    LT.derive$ex.derived <- LT.derive$Tx/LT.derive$lx
    LT.derive$ex.original <- select.country$ex
    LT.derive$diff <- LT.derive$ex.original-LT.derive$ex.derived
    LT.derive$Country <- country.select
    LT.derive$Edu <- edu.select
    LT.derive$Sex <- sex.select

    return(LT.derive[,c("Country","Edu","Sex","Age","px","lx","dx","Lx",
                        "Tx","ex.derived","ex.original","diff")])
}
```
The following code applies the function to all 16 European countries by educational attainment, stratified by sex.
```r
#these are the country codes
edu.countries <- c("BG","DK","EE","EL","HR","IT","HU", #CZ is currently not available
                   "PL","PT","RO","SI","SK","FI","SE","NO")


###Females###
out.females <- c()

for (country.select in edu.countries) {

    for (edu.select in c("higher","middle","lower")) {

        out.females <- rbind(out.females,my.function(country.select, edu.select, "F"))
}
}


###Males###
out.males <- c()

for (country.select in edu.countries) {

    for (edu.select in c("higher","middle","lower")) {

        out.males <- rbind(out.males,my.function(country.select, edu.select, "M"))
}
}
```
Finally, I plot the difference between the original $e_x$ and the derived $e_x$.
```r
par(mfrow=c(3,3))
for (edu in c("higher","middle","lower")) {
    plot(1,1, type="n", xlim=c(1,16), ylim=c(-0.2,0.2),
         main=paste("Females",edu,sep=" "), xlab="Countries",
         ylab="LE 30 original - LE30 derived")
    points(1:15,out.females$diff[out.females$Edu==edu & out.females$Age==30])
    text(1:15,out.females$diff[out.females$Edu==edu & out.females$Age==30], 1:16,
         label=out.females$Country[out.females$Edu==edu & out.females$Age==30])
}

for (edu in c("higher","middle","lower")) {
    plot(1,1, type="n", xlim=c(1,16), ylim=c(-0.2,0.2),
         main=paste("Males",edu,sep=" "), xlab="Countries",
         ylab="LE 30 original - LE30 derived")
    points(1:15,out.males$diff[out.males$Edu==edu & out.males$Age==30])
    text(1:15,out.males$diff[out.males$Edu==edu & out.males$Age==30], 1:16,
         label=out.males$Country[out.males$Edu==edu & out.males$Age==30])
}
```


## Complete life tables by age and education (stratified by women and men)
This prints all the age- and education-specific life tables (the output it omitted).

```r
library(knitr)

table.fun <- function(country.select) {
    
    print(
        kable(filter(out.females, Country==country.select & Edu=="higher"),
              digits=4, caption=paste("Life table for high-educated women in",
                                      country.select,", 2016",sep=" ")) 
        )
    print(
        kable(filter(out.females, Country==country.select & Edu=="middle"),
              digits=4, caption=paste("Life table for middle-educated women in",
                                      country.select,", 2016",sep=" ")) 
          )
    print(
        kable(filter(out.females, Country==country.select & Edu=="lower"),
              digits=4, caption=paste("Life table for low-educated women in",
                                      country.select,", 2016",sep=" "))
        )

    print(
        kable(filter(out.males, Country==country.select & Edu=="higher"),
              digits=4, caption=paste("Life table for high-educated men in",
                                      country.select,", 2016",sep=" ")) 
            )
    
    print(
        kable(filter(out.males, Country==country.select & Edu=="middle"),
              digits=4, caption=paste("Life table for middle-educated men in",
                                      country.select,", 2016",sep=" ")) 
            )
    print(
        kable(filter(out.males, Country==country.select & Edu=="lower"),
              digits=4, caption=paste("Life table for low-educated men in",
                                      country.select,", 2016",sep=" ")) 
            )    
}

for (country in edu.countries) {
    table.fun(country)  
}
```
## References
- Sauerberg, M. (2021). The imapact of population's educational attainment on Healthy Life Years in Europe: An empirical illustration of 16 European countries. SSM - Population Health, 15(100857).