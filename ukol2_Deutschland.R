library(tidyquant)
library(ggplot2)
library(ggfortify)
library(dynlm)
library(patchwork)
library(tstools)
library(stargazer)
library(strucchange)

# Stažení dat  
gdp_raw <- tq_get("CLVMNACSCAB1GQDE",
                  get  = "economic.data",
                  from = "1995-01-01",
                  to   = "2018-12-31")

u_raw <- tq_get("LRUNTTTTDEQ156S",
                get  = "economic.data",
                from = "1995-01-01",
                to   = "2018-12-31")


#--------------------------------
# 1.1
#  Převod na ts 
lGNP <- ts(log(gdp_raw$price), start = c(1995, 1), frequency = 4)
U   <- ts(u_raw$price,   start = c(1995, 1), frequency = 4)

#  Vizualizace původních řad 
p1 <- autoplot(lGNP, size = 1, colour = "red", facets = FALSE) +
  ggtitle("log - Reálný HDP Německa") +
  ylab("log(mil. EUR)") + xlab("")

p2 <- autoplot(U, size = 1, colour = "blue", facets = FALSE) +
  ggtitle("Míra nezaměstnanosti Německa") +
  ylab("%") + xlab("")

p1 / p2

Time <- ts(seq_len(length(lGNP)), start = c(1995, 1), frequency = 4)


# Endogenní identifikace zlomů
bp_lGNP <- breakpoints(lGNP ~ Time)
bp_U    <- breakpoints(U ~ Time)

summary(bp_lGNP)  # BIC optimum → 2008Q4
summary(bp_U)     # BIC monotónně klesá → volíme 2008Q4

par(mfrow = c(1, 2))
plot(bp_lGNP, main = "Zlomy — log HDP")
plot(bp_U,    main = "Zlomy — Nezaměstnanost")
par(mfrow = c(1, 1))


# Chowův test pro 2008Q4
break_point <- which(time(lGNP) == (2008 + 3/4))  # vrátí 56
sctest(lGNP ~ Time, type = "Chow", point = break_point)
sctest(U ~ Time,    type = "Chow", point = break_point)
#--------------------------------------------------------------
# 1.2 Tvorba umělých proměných a detrendování


# Dummy: 0 před 2008Q4, 1 od 2008Q4 dále
D <- create_dummy_ts(
  start_basic = c(1995, 1),
  end_basic   = c(2018, 4),
  dummy_start = c(2008, 4),
  dummy_end   = NULL,
  sp          = FALSE,
  basic_value = 0,
  dummy_value = 1,
  frequency   = 4
)

break_val <- window(Time, start = c(2008, 4), end = c(2008, 4))
Time_2    <- Time - break_val[1]

pom_data  <- ts.union(lGNP, U, Time, D, Time_2)

#  Regrese trendu pro log HDP (rovnice 5) 
trend_lGNP <- dynlm(lGNP ~ Time + D + I(D * Time_2),
                    data = pom_data)
summary(trend_lGNP)

#  Regrese trendu pro míru nezaměstnanosti (rovnice 6) 
trend_U <- dynlm(U ~ Time + D + I(D * Time_2),
                 data = pom_data)
summary(trend_U)

#---------------------------------------
# 1.3 Cyklické složky = residua z trendových regresí
y_cyc <- ts(residuals(trend_lGNP), start = c(1995, 1), frequency = 4)  # mezera výstupu
u_cyc <- ts(residuals(trend_U),    start = c(1995, 1), frequency = 4)  # mezera nezaměstnanosti

p3 <- autoplot(y_cyc, size = 1, colour = "red") +
  ggtitle("Mezera výstupu") +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype = "dashed")

p4 <- autoplot(u_cyc, size = 1, colour = "blue") +
  ggtitle("Mezera nezaměstnanosti") +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype = "dashed")

p3 / p4



#Vysledky
stargazer(trend_lGNP, trend_U, type = "text",
          title = "Trendové regrese — Německo 1995Q1–2018Q4",
          dep.var.labels = c("log HDP", "Míra nezaměstnanosti"),
          covariate.labels = c("Konstanta", "Čas", "Dummy 2008Q4", "Dummy × Čas od zlomu"),
          intercept.bottom = FALSE)


