#install.packages("ggthemes")
#install.packages("quantmod")
#install.packages("tseries")
#install.packages("ROI")
#install.packages("ROI.plugin.alabama")
#install.packages("PerformanceAnalytics")
#install.packages("ggridges")
#install.packages("xts")


#rm(list = ls(all.names = TRUE)) 


#Library
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggridges)
library(quantmod)
library(gtools)
library(tseries)
library(PerformanceAnalytics)
library(ROI)
library(alabama)
library(xts)

#Define Variables

  #Financial Industry
ticker_list<-c("AIG", "BAC", "C", "JPM", "AON","AFL","STT","HBAN","AXP")

      #c("AFL","STT","AXP","HBAN")
      #c("BAC","AIG","AON","HBAN")
      #c("HBAN","C","AIG","JPM")

  #  Health Care Industry
ticker_list<-c("DHR", "HUM", "JNJ", "LLY", "MDT","MRK","MYL","TFX","UNH")

      #c("UNH","TFX","DHR","MYL")
      #c("HUM","UNH","JNJ","DHR")
      #c("UNH","MRK","MDT","TFX")

  #Mix Industry
ticker_list<-c("WMB", "APA", "PG", "AAPL", "KO","HSY","LUV","PFE","GE")

        #c("APA","HSY","LUV","PG")
        #c("PG","PFE","AAPL","LUV")
        #c("KO","GE","AAPL","WMB")

#tickers <-sample(ticker_list,4)

tickers <-c("APA","HSY","LUV","PG")


tickers
fromdate<-"1985-01-01"
todate<-"1995-12-31"

fromdateEVA<-"1996-01-01"
todateEVA<-"2015-12-31"

#Data Download
n<-length(tickers)

portfolioPrices <- NULL 

for(t in tickers){
  
  portfolioPrices <-cbind(portfolioPrices, getSymbols(t,from=fromdate,to=todate,periodicity="monthly", auto.assign = FALSE)[,4])

  }

colnames(portfolioPrices)<-tickers

Portafolio<-ROC(portfolioPrices, type="discrete")[-1,]
head(Portafolio)


#Covariance Matrix calculations

covmatrix=matrix(c(cov(Portafolio)),nrow=n,ncol=n)
dimnames(covmatrix) = list(tickers, tickers)
covmatrix

#Correlation Matrix
corrmatrix<-cov2cor(covmatrix)
corrmatrix

chart.Correlation(Portafolio)

#Variables for PK

d=sqrt(diag(covmatrix))
D=diag(d,nrow=ncol(covmatrix),ncol=ncol(covmatrix))
Dinv<-solve(D)

#Dinv<-diag(1,nrow(D))
#  for(i in 1:nrow(D)){
#   Dinv[i,i]<-1/D[i,i]
#  }

PCA=eigen(covmatrix)
Lambda=sqrt(diag(PCA$values))
P=PCA$vectors


USV=svd(Lambda%*%t(P)%*%D)
USV
U=USV$u
S=diag(USV$d)
V=USV$v


#Function ENB for the weight Wk of the stocks in the portafolio


ENB<-function(w){
  
Ainvw<-Dinv%*%V%*%t(U)%*%Lambda%*%t(P)%*%w
WEW<-as.numeric(w%*%covmatrix%*%w)
    
Pk<-rep(1,n)
for(i in 1:n){
    
    Pk[i]<-(covmatrix[i,i]*Ainvw[i]^2)/WEW
    
  }
  
lnpk<-log(Pk)
pklnpk<-as.numeric(Pk%*%lnpk)
efnb<-exp(-1*pklnpk)
  
return(efnb)
  
}


#Optimization Model

sem<-rep(1/n,n)


Rest_i<-diag(rep(1,n+1),n+1)
Rest_i[1,]<-rep(1,n+1)
Rest_i<-Rest_i[,-1]
Rest_d

Rest_d<-c(1, rep(0,n))
Rest_d<-matrix(Rest_d,ncol=1)

Rest_v<-c("==",rep(">=",n))
Rest_v<-matrix(Rest_v,ncol=1)

constraints<-L_constraint(Rest_i,Rest_v,rhs=Rest_d)
Objective<-F_objective(ENB,n)
ProbOpt<-OP(Objective,constraints,maximum=TRUE)
Solve<-ROI_solve(ProbOpt,control=list(start=sem))
W_ENB<-solution(Solve)

format(W_ENB,scientific = F)
W_ENB

W_ENB<-round(W_ENB,5)
names(W_ENB) <- tickers

#Markowitz and NAIV

opt<- portfolio.optim(as.matrix(Portafolio))
W_mkwtz<-opt$pw
names(W_mkwtz) <- tickers


W_NAIV<-rep(1/n,n)
names(W_NAIV) <- tickers


#Weights Summary

weights<-cbind(as.data.frame(W_ENB),as.data.frame(W_mkwtz),as.data.frame(W_NAIV))
weights<-cbind(Stocks=rownames(weights),weights)
rownames(weights)<-1:nrow(weights)

wp<-gather(weights,"method","weight",2:4)

ggplot(wp,aes(Stocks, weight, color=method))+
  geom_point(size=5,alpha=0.7) +
  labs(title="Peso por Metodología",x="Acciones",y="Pesos",color="Metodología")+
  scale_color_tableau(labels = c("ENB", "Markowitz","NAIV"))+
  ylim(0,NA)+
  theme_gdocs()


#Data Evaluation

portfolioPricesEVA <- NULL 

for(t in tickers){
  
  portfolioPricesEVA <-cbind(portfolioPricesEVA, getSymbols(t,from=fromdateEVA,to=todateEVA,periodicity="monthly", auto.assign = FALSE)[,4])
  
}

colnames(portfolioPricesEVA)<-tickers

PortafolioEVA<-ROC(portfolioPricesEVA, type="discrete")[-1,]
head(portfolioPricesEVA)
head(PortafolioEVA)


ENB<-Return.portfolio(PortafolioEVA, weights=W_ENB,rebalance_on="years")
names(ENB)<-"ENB"
MKWTZ<-Return.portfolio(PortafolioEVA, weights=W_mkwtz,rebalance_on="years")
names(MKWTZ)<-"MKWTZ"
NAIV<-Return.portfolio(PortafolioEVA, weights=W_NAIV,rebalance_on="years")
names(NAIV)<-"NAIV"

# Returns

Methods<-cbind(ENB,MKWTZ,NAIV)
Methods_df<-as.data.frame(Methods)
Methods_df<-cbind(Fecha=rownames(Methods_df),Methods_df)
rownames(Methods_df)<-1:nrow(Methods_df)
head(Methods_df)

Methods_df$Fecha<-as.Date(Methods_df$Fecha)

Methods_plot<-gather(Methods_df,"Method","Return",2:4)
Methods_plot<-Methods_plot %>% mutate(Periodos=case_when(Fecha<"2001-12-31" ~ "Periodo 1996-2001",
                                                         Fecha >="2001-12-31" & Fecha<"2009-12-31"~ "Periodo 2002-2009",
                                                         Fecha >="2009-12-31"~ "Periodo 2010-2015"))
glimpse(Methods_plot)


    #Plot of Returns distribution per methodology
ggplot(Methods_plot,aes(Method,Return, color=Method))+
  geom_hline(yintercept=0,color="black")+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.2)+
  facet_grid(~Periodos)+
  labs(x="", y="Retornos Mensuales", title = "Distribución de retornos")+
  theme_gdocs()+
  scale_color_tableau(name="Metodología")

 
    #Wealth Returns

wealth_ENB<-Return.portfolio(PortafolioEVA, weights=W_ENB,rebalance_on="years",wealth.index = TRUE)
names(wealth_ENB)<-"ENB"
wealth_MKWTZ<-Return.portfolio(PortafolioEVA, weights=W_mkwtz,rebalance_on="years",wealth.index = TRUE)
names(wealth_MKWTZ)<-"MKWTZ"
wealth_NAIV<-Return.portfolio(PortafolioEVA, weights=W_NAIV,rebalance_on="years",wealth.index = TRUE)
names(wealth_NAIV)<-"NAIV"

wealth_Methods<-cbind(wealth_ENB,wealth_MKWTZ,wealth_NAIV)

wealth_Methods_df<-as.data.frame(wealth_Methods)
wealth_Methods_df<-cbind(Fecha=rownames(wealth_Methods_df),wealth_Methods_df)
rownames(wealth_Methods_df)<-1:nrow(wealth_Methods_df)
head(wealth_Methods_df)

wealth_Methods_df$Fecha<-as.Date(wealth_Methods_df$Fecha)
wealth_Methods_plot<-gather(wealth_Methods_df,"Method","Wealth",2:4)
glimpse(wealth_Methods_plot)

    #Plot of Wealth per methodology
ggplot(wealth_Methods_plot,aes(Fecha,Wealth,color=Method))+
  geom_line(size=1,alpha=0.7)+
  labs(title = "Inversión en el tiempo de 1$ por metodología")+
  theme_gdocs()+
  scale_color_tableau(name="Metodología")
 
#SharpeRatio

head(Methods_df)


Methods_sharpe<-Methods_df %>% group_by(Year=year(Fecha)) %>% 
  summarise(SharpeRatio_ENB=mean(ENB)/sd(ENB),
            SharpeRatio_MKWTZ=mean(MKWTZ)/sd(MKWTZ),
            SharpeRatio_NAIV=mean(NAIV)/sd(NAIV)) 

Methods_sharpe<-Methods_sharpe %>%  mutate(Winner= case_when(SharpeRatio_ENB>SharpeRatio_MKWTZ & SharpeRatio_ENB> SharpeRatio_NAIV~"ENB",
                           SharpeRatio_MKWTZ>SharpeRatio_ENB & SharpeRatio_MKWTZ>SharpeRatio_NAIV~"MKWTZ",
                           TRUE~"NAIV"
                            ))


table(Methods_sharpe$Winner)

    #Sharpie Ratio Plot

Methods_sharpe_plot<-gather(Methods_sharpe,"Method","Sharpe_Ratio",2:4)

tail(Methods_sharpe_plot)

ggplot(Methods_sharpe_plot,aes(as.factor(Year),Sharpe_Ratio,fill=Winner))+
  geom_col(alpha=0.3,size=4,position="dodge")+
  geom_point(aes(color=Method),size=2)+
  geom_hline(yintercept=0,color="black")+
  scale_fill_tableau(name="Ganador")+
  scale_color_tableau(name="Metodología", labels=c("ENB","MKWTZ","NAIV"))+
  labs(x="", y="Sharpe Ratio", title="Desempeño y ganador según Sharpe Ratio")+
  theme_gdocs()
  

#Standard Deviation

Methods_sd<-Methods_df %>% group_by(Year=year(Fecha)) %>% 
  summarise(sd_ENB=sd(ENB),
            sd_MKWTZ=sd(MKWTZ),
            sd_NAIV=sd(NAIV)) 

Methods_sd<-Methods_sd %>%  mutate(Winner= case_when(sd_ENB<sd_MKWTZ & sd_ENB< sd_NAIV~"ENB",
                                                             sd_MKWTZ<sd_ENB & sd_MKWTZ<sd_NAIV~"MKWTZ",
                                                             TRUE~"NAIV"))

table(Methods_sd$Winner)

    #Standar Deviation Plot 

Methods_sd_plot<-gather(Methods_sd,"Method","sd",2:4)

ggplot(Methods_sd_plot, aes(x=as.factor(Year),y=sd,fill=Winner))+
  geom_col(alpha=0.3,size=4,position="dodge")+
  geom_point(aes(color=Method),size=2)+
  scale_fill_tableau(name="Ganador")+
  scale_color_tableau(name="Metodología", labels=c("ENB","MKWTZ","NAIV"))+
  labs(x="", y="Desviación Estandar", title="Desempeño y ganador según Desviación Estandar")+
  theme_gdocs()
 

#Turnover

    #ENB turnover

ENB_TO<-Return.portfolio(PortafolioEVA,rebalance_on="years", weights=W_ENB,verbose=TRUE)

beginWeights <- ENB_TO$BOP.Weight
endWeights <- ENB_TO$EOP.Weight
txns <- beginWeights - lag.xts(endWeights)
monthlyTO_ENB <- xts(rowSums(abs(txns[,1:4])), order.by=index(txns))


    #MKWTZ turnover

MKWTZ_TO<-Return.portfolio(PortafolioEVA, weights=W_mkwtz,rebalance_on="years",verbose=TRUE)

beginWeights2 <- MKWTZ_TO$BOP.Weight
endWeights2 <- MKWTZ_TO$EOP.Weight
txns2 <- beginWeights2 - lag.xts(endWeights2)
monthlyTO_MKWTZ <- xts(rowSums(abs(txns2[,1:4])), order.by=index(txns2))

    #NAIV turnover

NAIV_TO<-Return.portfolio(PortafolioEVA, weights=W_NAIV,rebalance_on="years",verbose=TRUE)

beginWeights3 <- NAIV_TO$BOP.Weight
endWeights3 <- NAIV_TO$EOP.Weight
txns3 <- beginWeights3 - lag.xts(endWeights3)
monthlyTO_NAIV <- xts(rowSums(abs(txns3[,1:4])), order.by=index(txns3))


head(beginWeights3,24)
head(endWeights3,24)
head(Methods_TO_df,42)

    #Turnover plot

Methods_turnover<-cbind(monthlyTO_ENB,monthlyTO_MKWTZ,monthlyTO_NAIV)
Methods_TO_df<-as.data.frame(Methods_turnover)
Methods_TO_df<-cbind(Fecha=rownames(Methods_TO_df),Methods_TO_df)
rownames(Methods_TO_df)<-1:nrow(Methods_TO_df)

Methods_TO_df$Fecha<-as.Date(Methods_df$Fecha)

Methods_TO_df<-Methods_TO_df %>%  mutate(Winner= case_when(monthlyTO_ENB<monthlyTO_MKWTZ & monthlyTO_ENB< monthlyTO_NAIV~"ENB",
                                                       monthlyTO_MKWTZ<monthlyTO_ENB & monthlyTO_MKWTZ<monthlyTO_NAIV~"MKWTZ",
                                                       monthlyTO_NAIV<monthlyTO_ENB & monthlyTO_NAIV<monthlyTO_MKWTZ~"NAIV",
                                                       TRUE~"NA"))

      #Turnover Plot 

Methods_TO_plot<-gather(Methods_TO_df,"Method","Turnover",2:4)
head(Methods_TO_plot,2)


Methods_TO_plot<-Methods_TO_plot %>% 
    group_by(as.factor(year(Fecha)),Method,Winner) %>% 
    summarize(sum(Turnover)) %>% 
    filter(Winner != "NA")


head(Methods_TO_plot,2)
colnames(Methods_TO_plot)<-c("year","Method","Winner","Turnover")
Methods_TO_plot<-as.data.frame(Methods_TO_plot)


table(Methods_TO_plot$Winner)/3

ggplot(Methods_TO_plot, aes(x=year,y=Turnover,fill=Winner))+
  geom_col(alpha=0.3,size=4,position="dodge")+
  geom_point(aes(color=Method),alpha=0.6,size=2)+
  scale_fill_tableau(name="Ganador")+
  scale_color_tableau(name="Metodología", labels=c("ENB","MKWTZ","NAIV"))+
  labs(x="",y="Turnover", title="Desempeño y ganador según Turnover")+
  theme_gdocs()
 


# Resumen de Tablas

  #Sharpie Ratio
  table(Methods_sharpe$Winner)

  #SD
  table(Methods_sd$Winner)

  #Turnover
  table(Methods_TO_plot$Winner)/3



