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

# PŘÍKLAD 2: Odhad Okunova vztahu

library(car)      
library(ggplot2)
library(patchwork)
library(stargazer)
library(strucchange)


# data

df2 <- data.frame(
  u_cyc = as.numeric(u_cyc),
  y_cyc = as.numeric(y_cyc),
  D     = as.numeric(D),
  row.names = NULL
)
df2$D_y   <- df2$D * df2$y_cyc          # interakce pro (3')
df2$y_pos <- pmax(df2$y_cyc, 0)         # kladná mezera výstupu
df2$y_neg <- pmin(df2$y_cyc, 0)         # záporná mezera výstupu

n_total <- nrow(df2)                     # 96 čtvrtletí
n1      <- break_point - 1              # obs. PŘED zlomem
n2      <- n_total - n1                 # obs. OD zlomu   

# 2.1  Základní statická verze Okunova vztahu – rovnice (3)
#      u_cyc,t = beta * y_cyc,t + eps_t

okun_basic <- lm(u_cyc ~ y_cyc, data = df2)
summary(okun_basic)

# 2.2  Testování strukturální stability koeficientů

# a) Standardní Chow F-test
chow_std <- sctest(u_cyc ~ y_cyc, type = "Chow",
                   point = break_point, data = df2)
cat("\n--- Standardní Chowův test (2008Q4) ---\n")
print(chow_std)

# b) forecast Chow test

okun_pre_fit <- lm(u_cyc ~ y_cyc, data = df2[seq_len(n1), ])

RSS1   <- sum(residuals(okun_pre_fit)^2)
u_fc   <- predict(okun_pre_fit,
                  newdata = df2[(n1 + 1):n_total, , drop = FALSE])
RSS_fc <- sum((df2$u_cyc[(n1 + 1):n_total] - u_fc)^2)
k_par  <- length(coef(okun_pre_fit))          # 2 (intercept + sklon)

F_fc <- (RSS_fc / n2) / (RSS1 / (n1 - k_par))
p_fc <- pf(F_fc, n2, n1 - k_par, lower.tail = FALSE)

cat(sprintf(
  "\n--- Chow forecast test ---\nF(%d, %d) = %.4f,  p-hodnota = %.4g\n",
  n2, n1 - k_par, F_fc, p_fc
))


# 2.3  Okunův vztah se strukturálním zlomem – rovnice (3')
#      u_cyc = beta0 + beta1*y_cyc + delta0*D + delta1*(D*y_cyc) + eps

okun_break <- lm(u_cyc ~ y_cyc + D + D_y, data = df2)
summary(okun_break)

# Okunovy koeficienty před a po zlomu
beta_pre  <- coef(okun_break)["y_cyc"]
beta_post <- coef(okun_break)["y_cyc"] + coef(okun_break)["D_y"]
cat(sprintf("\nOkunův koeficient PŘED zlomem : %.4f\n", beta_pre))
cat(sprintf("Okunův koeficient PO  zlomu   : %.4f\n", beta_post))

# Test, zda se sklon skutečně změnil (H0: delta1 = 0)
cat("\nTest změny sklonu (delta1 = 0):\n")
print(summary(okun_break)$coefficients["D_y", ])

# 2.4  Asymetrie Okunova koeficientu
#      u_cyc = beta_pos * y_pos + beta_neg * y_neg + eps
okun_asym <- lm(u_cyc ~ y_pos + y_neg, data = df2)
summary(okun_asym)

# waldův test H0: beta_+ = beta_-
asym_test <- linearHypothesis(okun_asym, "y_pos = y_neg")
cat("\n--- Waldův test asymetrie (beta_+ = beta_-) ---\n")
print(asym_test)

# vizualizace
# scatter
p_basic <- ggplot(df2, aes(x = y_cyc, y = u_cyc)) +
  geom_point(colour = "steelblue", alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "firebrick", se = TRUE, linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "Okunův zákon – základní verze",
       x = "Mezera výstupu",
       y = "Mezera nezaměstnanosti") +
  theme_minimal(base_size = 11)

# scatter: před a po zlomu
df2$period <- factor(
  ifelse(df2$D == 0, "1995Q1\u20132008Q3", "2008Q4\u20132018Q4"),
  levels = c("1995Q1\u20132008Q3", "2008Q4\u20132018Q4")
)

p_break <- ggplot(df2, aes(x = y_cyc, y = u_cyc,
                           colour = period, shape = period)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  labs(title = "Okunův zákon – strukturální zlom",
       x = "Mezera výstupu",
       y = "Mezera nezaměstnanosti",
       colour = "Období", shape = "Období") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_basic / p_break   # vykreslení přes patchwork


# tabulky výsledků
stargazer(okun_basic, okun_break, okun_asym,
          type      = "text",
          title     = "Odhady Okunova vztahu — Německo 1995Q1–2018Q4",
          dep.var.labels   = "Mezera nezaměstnanosti",
          column.labels    = c("Základ. (3)", "Se zlomem (3')", "Asymetr."),
          covariate.labels = c("Mezera výstupu (\\hat{y})",
                               "Dummy 2008Q4 (D)",
                               "$D \\times \\hat{y}$",
                               "Mezera výstupu kladná",
                               "Mezera výstupu záporná",
                               "Konstanta"),
          intercept.bottom = TRUE,
          omit.stat = c("ser", "f"),
          notes     = "Směrodatné chyby v závorkách. *p<0.1; **p<0.05; ***p<0.01")

# PŘÍKLAD 3: Dynamický Okunův koeficient

library(car)
library(MASS)      # mvrnorm pro simulační metodu
library(stargazer)


# Příprava dat — tvorba posunutých (lagged) proměnných

n      <- length(u_cyc)
lag_k  <- function(x, k) c(rep(NA, k), head(as.numeric(x), n - k))

MAX_P  <- 4   # maximální AR řád
MAX_Q  <- 4   # maximální zpoždění mezery výstupu

df3 <- data.frame(
  u = as.numeric(u_cyc),
  y = as.numeric(y_cyc),
  D = as.numeric(D)
)

# Zpožděné proměnné: uLk = L^k(u), yLk = L^k(y)
for (k in seq_len(MAX_P)) df3[[paste0("uL", k)]]  <- lag_k(df3$u, k)
for (k in seq(0, MAX_Q))  df3[[paste0("yL", k)]]  <- lag_k(df3$y, k)

# Interakce D × L^k(y)  a  D × L^k(u) — pro modely se zlomem
for (k in seq(0, MAX_Q))  df3[[paste0("DyL", k)]] <- df3$D * lag_k(df3$y, k)
for (k in seq_len(MAX_P)) df3[[paste0("DuL", k)]] <- df3$D * lag_k(df3$u, k)


# 3.1  Výběr optimálního počtu zpoždění — rovnice (7)
#  û_t = α_0 + Σ_{i=1}^p α_i û_{t-i} + Σ_{j=0}^q β_j ŷ_{t-j} + ε_t


ic_grid <- data.frame(p=integer(), q=integer(), AIC=numeric(), BIC=numeric())

for (p in seq_len(MAX_P)) {
  for (q in seq(0, MAX_Q)) {
    rhs <- paste(c(paste0("uL", seq_len(p)),
                   paste0("yL", seq(0, q))),
                 collapse = " + ")
    m <- tryCatch(
      lm(as.formula(paste("u ~", rhs)), data = df3, na.action = na.omit),
      error = function(e) NULL
    )
    if (!is.null(m))
      ic_grid <- rbind(ic_grid,
                       data.frame(p=p, q=q, AIC=AIC(m), BIC=BIC(m)))
  }
}

best_bic <- ic_grid[which.min(ic_grid$BIC), ]
best_aic <- ic_grid[which.min(ic_grid$AIC), ]
cat("BIC optimum: p =", best_bic$p, ", q =", best_bic$q, "\n")
cat("AIC optimum: p =", best_aic$p, ", q =", best_aic$q, "\n")

# Pomocná funkce pro odhad modelu (7) s volitelným rozšířením pravé strany
fit_dyn <- function(p, q, extra = NULL) {
  rhs <- paste(c(paste0("uL", seq_len(p)),
                 paste0("yL", seq(0, q)),
                 extra),
               collapse = " + ")
  lm(as.formula(paste("u ~", rhs)), data = df3, na.action = na.omit)
}

# Optimální modely dle BIC a AIC
p_opt <- best_bic$p; q_opt <- best_bic$q

okun_dyn     <- fit_dyn(p_opt, q_opt)                    
okun_dyn_aic <- fit_dyn(best_aic$p, best_aic$q)          

summary(okun_dyn)
summary(okun_dyn_aic)

# Výpočet dlouhodobého Okunova koeficientu
lr_coef <- function(model, p, q,
                    break_alpha = FALSE, break_beta = FALSE) {
  cf <- coef(model)
  sa <- sum(cf[paste0("uL", seq_len(p))])
  sb <- sum(cf[paste0("yL", seq(0, q))])

  if (break_alpha) {
    du_nms <- paste0("DuL", seq_len(p))
    sa <- sa + sum(cf[du_nms[du_nms %in% names(cf)]])
  }
  if (break_beta) {
    dy_nms <- paste0("DyL", seq(0, q))
    sb <- sb + sum(cf[dy_nms[dy_nms %in% names(cf)]])
  }
  sb / (1 - sa)
}

lr7_bic <- lr_coef(okun_dyn,     p_opt,        q_opt)
lr7_aic <- lr_coef(okun_dyn_aic, best_aic$p,   best_aic$q)
cat(sprintf("LR (BIC model, p=%d q=%d): %.3f\n", p_opt, q_opt, lr7_bic))
cat(sprintf("LR (AIC model, p=%d q=%d): %.3f\n", best_aic$p, best_aic$q, lr7_aic))

# 3.2  Strukturální zlom v rovnici (7') — varianty a) a b)
#  Varianta a): zlom pouze v β  (koeficienty u ŷ_{t-j})
#    û_t = α_0 + Σ α_i û_{t-i} + Σ (β_j + δ_j D_t) ŷ_{t-j} + ε_t
#
#  Varianta b): zlom v α i β
#    û_t = α_0 + Σ (α_i + γ_i D_t) û_{t-i} + Σ (β_j + δ_j D_t) ŷ_{t-j} + ε_t

okun_7a <- fit_dyn(p_opt, q_opt,
                   extra = paste0("DyL", seq(0, q_opt)))

okun_7b <- fit_dyn(p_opt, q_opt,
                   extra = c(paste0("DyL", seq(0, q_opt)),
                              paste0("DuL", seq_len(p_opt))))

summary(okun_7a)
summary(okun_7b)

# Chowovy testy strukturální stability
chow_7a <- linearHypothesis(okun_7a,
             paste0(paste0("DyL", seq(0, q_opt)), " = 0"))
chow_7b <- linearHypothesis(okun_7b,
             c(paste0(paste0("DyL", seq(0, q_opt)), " = 0"),
               paste0(paste0("DuL", seq_len(p_opt)), " = 0")))

cat("\n--- Chowův test — varianta a) (zlom v β) ---\n")
print(chow_7a)
cat("\n--- Chowův test — varianta b) (zlom v α a β) ---\n")
print(chow_7b)

# Dlouhodobé koeficienty před a po zlomu
lr_7a_pre  <- lr_coef(okun_7a, p_opt, q_opt, FALSE, FALSE)
lr_7a_post <- lr_coef(okun_7a, p_opt, q_opt, FALSE, TRUE)
lr_7b_pre  <- lr_coef(okun_7b, p_opt, q_opt, FALSE, FALSE)
lr_7b_post <- lr_coef(okun_7b, p_opt, q_opt, TRUE,  TRUE)

cat(sprintf("\nVar. a)  LR před: %.3f  |  LR po: %.3f\n", lr_7a_pre, lr_7a_post))
cat(sprintf("Var. b)  LR před: %.3f  |  LR po: %.3f\n", lr_7b_pre, lr_7b_post))

# ---------------------------------------------------------------
# 3.3  Plný strukturální zlom — všechny koeficienty včetně konstanty
# ---------------------------------------------------------------

okun_full <- fit_dyn(p_opt, q_opt,
                     extra = c("D",
                               paste0("DyL", seq(0, q_opt)),
                               paste0("DuL", seq_len(p_opt))))
summary(okun_full)

chow_full <- linearHypothesis(okun_full,
               c("D = 0",
                 paste0(paste0("DyL", seq(0, q_opt)), " = 0"),
                 paste0(paste0("DuL", seq_len(p_opt)), " = 0")))
cat("\n--- Chowův test — plný zlom ---\n")
print(chow_full)

lr_full_pre  <- lr_coef(okun_full, p_opt, q_opt, FALSE, FALSE)
lr_full_post <- lr_coef(okun_full, p_opt, q_opt, TRUE,  TRUE)
cat(sprintf("Plný zlom  LR před: %.3f  |  LR po: %.3f\n", lr_full_pre, lr_full_post))


# 3.4  Směrodatné odchylky LR koeficientů
# Deltova metoda a simulační (parametrický bootstrap)


B <- 10000   # počet Monte Carlo výběrů
set.seed(42)

# Numerický gradient pro Deltovu metodu 
delta_se <- function(model, lr_fn, h = 1e-6) {
  cf   <- coef(model)
  Sig  <- vcov(model)
  lr0  <- lr_fn(cf)
  grad <- vapply(seq_along(cf), function(i) {
    cf2    <- cf
    cf2[i] <- cf2[i] + h
    (lr_fn(cf2) - lr0) / h
  }, numeric(1))
  sqrt(as.numeric(t(grad) %*% Sig %*% grad))
}

# --- Simulační metoda ---
sim_se <- function(model, lr_fn, B = 10000) {
  cf    <- coef(model)
  Sig   <- vcov(model)
  draws <- MASS::mvrnorm(B, mu = cf, Sigma = Sig)
  lrs   <- apply(draws, 1, function(th) {
    names(th) <- names(cf)
    tryCatch(lr_fn(th), error = function(e) NA_real_)
  })
  sd(lrs, na.rm = TRUE)
}

# Uzávěry (closures) pro každý LR koeficient 
make_lr <- function(p, q, ba = FALSE, bb = FALSE) {
  function(theta) {
    sa <- sum(theta[paste0("uL", seq_len(p))])
    sb <- sum(theta[paste0("yL", seq(0, q))])
    if (ba) {
      du_nms <- paste0("DuL", seq_len(p))
      sa <- sa + sum(theta[du_nms[du_nms %in% names(theta)]])
    }
    if (bb) {
      dy_nms <- paste0("DyL", seq(0, q))
      sb <- sb + sum(theta[dy_nms[dy_nms %in% names(theta)]])
    }
    sb / (1 - sa)
  }
}

# Seznam modelů a příslušných LR funkcí
specs <- list(
  list(model=okun_dyn,  fn=make_lr(p_opt,q_opt,F,F),  label="Základní (BIC)"),
  list(model=okun_7a,   fn=make_lr(p_opt,q_opt,F,F),  label="Zlom-β, před"),
  list(model=okun_7a,   fn=make_lr(p_opt,q_opt,F,T),  label="Zlom-β, po"),
  list(model=okun_7b,   fn=make_lr(p_opt,q_opt,F,F),  label="Zlom-αβ, před"),
  list(model=okun_7b,   fn=make_lr(p_opt,q_opt,T,T),  label="Zlom-αβ, po"),
  list(model=okun_full, fn=make_lr(p_opt,q_opt,F,F),  label="Plný zlom, před"),
  list(model=okun_full, fn=make_lr(p_opt,q_opt,T,T),  label="Plný zlom, po")
)

se_tbl <- do.call(rbind, lapply(specs, function(s) {
  lr_val <- s$fn(coef(s$model))
  data.frame(
    Model    = s$label,
    LR       = round(lr_val,  3),
    SE_Delta = round(delta_se(s$model, s$fn), 3),
    SE_Sim   = round(sim_se(s$model,   s$fn, B), 3),
    stringsAsFactors = FALSE
  )
}))

cat("\n=== Směrodatné odchylky LR koeficientů ===\n")
print(se_tbl, row.names = FALSE)

# Stargazer výstup základní tabulky modelů 
stargazer(okun_dyn, okun_7a, okun_7b, okun_full,
          type           = "text",
          title          = "Dynamické modely Okunova vztahu — Německo 1995Q1–2018Q4",
          dep.var.labels = "Mezera nezaměstnanosti",
          column.labels  = c("Základ. (7)", "Zlom-β (7'a)",
                             "Zlom-αβ (7'b)", "Plný zlom"),
          omit.stat      = c("ser", "f"),
          notes          = paste0("Sm. chyby v závorkách. ",
                                  "*p<0.1; **p<0.05; ***p<0.01"))
