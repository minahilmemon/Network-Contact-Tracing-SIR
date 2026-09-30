# ==============================================================
# NONLINEAR AGE-STRUCTURED NETWORK PDE
# WITH CONTACT TRACING
#
# GOLD-STANDARD probability equations; projection factors removed
#
# R_kj(t,a) = P_{I_j}^{S_k}(t,a) / Lambda_k(t)
#
# This avoids directly evaluating the small/small boundary
# P_{I_j}^{S_k}/Lambda_k.
#
# ABM is used only afterward as an independent check.
# ==============================================================

rm(list = ls())


# ==============================================================
# 1. PARAMETERS
# ==============================================================

N <- 100000

mean_degree <- 5

beta  <- 0.10
gamma <- 0.1

# tau = diagnosis rate
alpha <- 0.15

# theta = tracing-completion / tracing rate
delta <- 1.0

p <- 0.3

delta_p <- delta * p

I0 <- 20


# ==============================================================
# 2. NUMERICAL GRID
# ==============================================================

# Characteristic grid:
# dt = da

dt <- 0.2
da <- 0.2

t_max <- 200
a_max <- 30

times <- seq(
  0,
  t_max,
  by = dt
)

ages <- seq(
  0,
  a_max,
  by = da
)

Nt <- length(times)
Na <- length(ages)


# ==============================================================
# COMPOSITE SIMPSON RULE FOR AGE INTEGRALS
# ==============================================================
# The age grid is uniform. Composite Simpson's rule requires an
# even number of subintervals. With a_max = 30 and da = 0.2,
# there are 150 subintervals, so the rule applies exactly.

simpson_age <- function(y, h = da) {
  n_intervals <- length(y) - 1

  if (n_intervals < 2) {
    stop("Simpson's rule requires at least two subintervals.")
  }

  if (n_intervals %% 2 != 0) {
    stop("Simpson's rule requires an even number of subintervals.")
  }

  odd_sum <- sum(y[seq(2, length(y) - 1, by = 2)])
  even_sum <- sum(y[seq(3, length(y) - 2, by = 2)])

  (h / 3) * (
    y[1] +
      y[length(y)] +
      4 * odd_sum +
      2 * even_sum
  )
}


# ==============================================================
# 3. DEGREE DISTRIBUTION
# ==============================================================

# Junling's ABM uses a Poisson contact degree
# distribution with mean 5.

Kmax <- 20

all_degrees <- 0:Kmax

pk_all <- dpois(
  all_degrees,
  lambda = mean_degree
)

# Renormalize the truncated Poisson distribution.

pk_all <- pk_all / sum(pk_all)

# Degree-zero group

p0 <- pk_all[1]

# Positive-degree groups

degrees <- 1:Kmax

pk <- pk_all[-1]

K <- length(degrees)


# Size-biased degree distribution for a randomly
# chosen network neighbour:
#
# q_j = j p_j / <k>

mean_k <- sum(
  degrees * pk
)

q <- degrees * pk / mean_k

q <- q / sum(q)


# ==============================================================
# 4. INITIAL EPIDEMIC STATE
# ==============================================================

epsilon <- I0 / N

# Susceptible population by degree

S <- N * pk * (1 - epsilon)


# --------------------------------------------------------------
# Infection-age density I_k(t,a)
# --------------------------------------------------------------

I_age <- matrix(
  0,
  nrow = K,
  ncol = Na
)

# Initially the infected seeds are placed at infection age 0.

I_age[, 1] <- I0 * pk / da

# Degree-zero initially infected individuals

I_zero_0 <- I0 * p0


# ==============================================================
# 5. LAMBDA_k
#
# Lambda_k(t)
# =
# integral sum_j P_{I_j}^{S_k}(t,a) da
# ==============================================================

# For random initial infection, every susceptible focal
# degree class initially sees the same infected fraction.

Lambda <- rep(
  epsilon,
  K
)


# ==============================================================
# 6. NORMALIZED INFECTED-NEIGHBOUR DISTRIBUTION
#
# R[k,j,a]
# =
# P_{I_j}^{S_k}(t,a) / Lambda_k(t)
#
# Therefore:
#
# sum_j integral R[k,j,a] da = 1
# ==============================================================

R <- array(
  0,
  dim = c(K, K, Na)
)

for (ik in 1:K) {
  
  for (jj in 1:K) {
    
    R[ik, jj, 1] <-
      q[jj] / da
    
  }
  
}


# ==============================================================
# 7. P_{S_j}^{S_k}
# ==============================================================

PSS <- matrix(
  0,
  nrow = K,
  ncol = K
)

for (ik in 1:K) {
  
  PSS[ik, ] <-
    q * (1 - epsilon)
  
}


# ==============================================================
# 8. P_{T_j}^{S_k}
# ==============================================================

PST <- matrix(
  0,
  nrow = K,
  ncol = K
)


# ==============================================================
# 9. P_{S_j}^{I_k}(t,a)
# ==============================================================

PIS <- array(
  0,
  dim = c(K, K, Na)
)

for (ik in 1:K) {
  
  for (jj in 1:K) {
    
    PIS[ik, jj, 1] <-
      q[jj] * (1 - epsilon)
    
  }
  
}


# ==============================================================
# 10. P_{T_j}^{I_k}(t,a)
# ==============================================================

PTI <- array(
  0,
  dim = c(K, K, Na)
)


# ==============================================================
# 11. P_{I_j}^{I_k}(t,a1,a2)
# ==============================================================

PII <- array(
  0,
  dim = c(
    K,
    K,
    Na,
    Na
  )
)

# Random initial infected neighbours.
#
# This is an INITIAL CONDITION, not the infection-age
# boundary condition.

for (ik in 1:K) {
  
  for (jj in 1:K) {
    
    PII[
      ik,
      jj,
      1,
      1
    ] <-
      q[jj] *
      epsilon /
      da
    
  }
  
}


# ==============================================================
# 12. STORAGE
# ==============================================================

I_total <- numeric(Nt)

# At t = 0 the infected seeds are represented as a point mass at
# age zero (I0 * pk / da), not as a smooth age density. Therefore
# initialize the total directly rather than applying Simpson's rule
# to that discrete point mass.
I_total[1] <-
  I0 * sum(pk) +
  I_zero_0


Lambda_history <- matrix(
  0,
  nrow = Nt,
  ncol = K
)

Lambda_history[1, ] <- Lambda




# ==============================================================
# 13. START TIME LOOP
# ==============================================================

cat("\n")
cat("============================================\n")
cat("Starting GOLD-STANDARD nonlinear PDE simulation (dt=da=0.2, Simpson age integration)\n")
cat("============================================\n\n")


for (n in 1:(Nt - 1)) {
  
  
  # ------------------------------------------------------------
  # Progress display
  # ------------------------------------------------------------
  
  if (n %% 10 == 1) {
    
    cat(
      "PDE time =",
      times[n],
      "of",
      t_max,
      "\n"
    )
    
  }
  
  
  # ============================================================
  # A. TOTAL TRACED-NEIGHBOUR PROBABILITY
  #
  # P_T^{I_k}(t,a)
  # =
  # sum_j P_{T_j}^{I_k}(t,a)
  # ============================================================
  
  PT_total <- matrix(
    0,
    nrow = K,
    ncol = Na
  )
  
  for (ik in 1:K) {
    
    for (aa in 1:Na) {
      
      PT_total[ik, aa] <-
        sum(
          PTI[
            ik,
            ,
            aa
          ]
        )
      
    }
    
  }
  
  # ============================================================
  # B. W_kj(t,a)
  #
  # W_kj
  # =
  # integral P_{I_j}^{I_k}(t,a,a2) da2
  # ============================================================
  
  W <- array(
    0,
    dim = c(K, K, Na)
  )
  
  for (ik in 1:K) {
    
    for (jj in 1:K) {
      
      temp <-
        PII[
          ik,
          jj,
          ,
          
        ]
      
      W[
        ik,
        jj,
        
      ] <-
        apply(
          temp,
          1,
          simpson_age
        )
      
    }
    
  }
  
  
  # ============================================================
  # C. SUSCEPTIBLE POPULATION
  #
  # S'_k
  # =
  # - beta k S_k Lambda_k
  # ============================================================
  
  incidence <- numeric(K)
  
  for (ik in 1:K) {
    
    k <- degrees[ik]
    
    incidence[ik] <-
      beta *
      k *
      S[ik] *
      Lambda[ik]
    
  }
  
  
  S_new <-
    S -
    dt * incidence
  
  S_new <- pmax(
    S_new,
    0
  )
  
  
  # ============================================================
  # D. INFECTED POPULATION PDE
  #
  # (d_t + d_a) I_k
  #
  # =
  #
  # -[
  #    alpha
  #    + gamma
  #    + delta p k P_T^{I_k}
  #  ] I_k
  # ============================================================
  
  I_new <- matrix(
    0,
    nrow = K,
    ncol = Na
  )
  
  
  for (ik in 1:K) {
    
    k <- degrees[ik]
    
    loss <-
      alpha +
      gamma +
      delta_p *
      k *
      PT_total[ik, ]
    
    
    I_new[
      ik,
      2:Na
    ] <-
      I_age[
        ik,
        1:(Na - 1)
      ] *
      exp(
        -loss[
          1:(Na - 1)
        ] *
          dt
      )
    
    
    # Infection-age-zero boundary:
    #
    # I_k(t,0)
    # =
    # beta k S_k Lambda_k
    
    I_new[
      ik,
      1
    ] <-
      incidence[ik]
    
  }
  
  
  # ============================================================
  # E. UPDATE P_{S_j}^{S_k}
  # ============================================================
  
  PSS_new <- PSS
  
  
  for (ik in 1:K) {
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      rhs <-
        -beta *
        (j - 1) *
        PSS[ik, jj] *
        Lambda[jj] +
        beta *
        PSS[ik, jj] *
        Lambda[ik]
      
      
      PSS_new[
        ik,
        jj
      ] <-
        PSS[
          ik,
          jj
        ] +
        dt * rhs
      
    }
    
  }
  
  
  PSS_new <- pmax(
    PSS_new,
    0
  )
  
  
  # ============================================================
  # F. UPDATE P_{T_j}^{S_k}
  #
  # Since
  #
  # P_{I_j}^{S_k}
  # =
  # Lambda_k R_kj,
  #
  # its integrals are calculated without forming a
  # tiny numerator explicitly.
  # ============================================================
  
  PST_new <- PST
  
  
  for (ik in 1:K) {
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      
      integral_R <-
        simpson_age(
          R[
            ik,
            jj,
            
          ]
        )
      
      
      tracing_average <-
        simpson_age(
          R[
            ik,
            jj,
            
          ] *
            PT_total[
              jj,
              
            ]
        )
      
      
      rhs <-
        alpha *
        Lambda[ik] *
        integral_R +
        delta_p *
        (j - 1) *
        Lambda[ik] *
        tracing_average +
        beta *
        PST[
          ik,
          jj
        ] *
        Lambda[ik]
      
      
      PST_new[
        ik,
        jj
      ] <-
        PST[
          ik,
          jj
        ] +
        dt * rhs
      
    }
    
  }
  
  
  PST_new <- pmax(
    PST_new,
    0
  )
  
  
  # ============================================================
  # G. EVOLVE Lambda_k
  #
  # Derived by integrating the
  # P_{I_j}^{S_k} transport equation over age
  # and summing over j.
  # ============================================================
  
  Lambda_dot <- numeric(K)
  
  
  for (ik in 1:K) {
    
    
    # Infection-age-zero influx into susceptible-infected
    # neighbour probabilities.
    
    boundary_influx <- 0
    
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      boundary_influx <-
        boundary_influx +
        beta *
        (j - 1) *
        PSS[
          ik,
          jj
        ] *
        Lambda[jj]
      
    }
    
    
    # Tracing-weighted average appearing in the
    # integrated transport equation.
    
    tracing_average <- 0
    
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      tracing_average <-
        tracing_average +
        (j - 1) *
        simpson_age(
          R[
            ik,
            jj,
            
          ] *
            PT_total[
              jj,
              
            ]
        )
      
    }
    
    
    Lambda_dot[ik] <-
      boundary_influx +
      Lambda[ik] *
      (
        -(beta + alpha + gamma) +
          beta * Lambda[ik] -
          delta_p *
          tracing_average
      )
    
  }
  
  
  Lambda_new <-
    Lambda +
    dt *
    Lambda_dot
  
  
  # We are solving the NONLINEAR epidemic away from DFE.
  #
  # This tiny floor is only a machine-arithmetic safeguard.
  # It is NOT a declaration that the DFE boundary equals zero.
  
  Lambda_new <- pmax(
    Lambda_new,
    1e-14
  )

  
  
  # ============================================================
  # I. EVOLVE NORMALIZED R_kj
  #
  # P_{I_j}^{S_k}
  # =
  # Lambda_k R_kj
  #
  # Therefore
  #
  # (d_t+d_a)R
  #
  # =
  #
  # [coefficient - Lambda'_k/Lambda_k]R
  #
  # This is where the common small epidemic amplitude has
  # analytically cancelled.
  # ============================================================
  
  R_new <- array(
    0,
    dim = c(K, K, Na)
  )
  
  
  for (ik in 1:K) {
    
    
    growth_normalization <-
      Lambda_dot[ik] /
      max(
        Lambda[ik],
        1e-14
      )
    
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      
      rate_R <-
        -(beta + alpha + gamma) -
        delta_p *
        (j - 1) *
        PT_total[
          jj,
          
        ] +
        beta *
        Lambda[ik] -
        growth_normalization
      
      
      # Characteristic transport
      
      R_new[
        ik,
        jj,
        2:Na
      ] <-
        R[
          ik,
          jj,
          1:(Na - 1)
        ] *
        exp(
          rate_R[
            1:(Na - 1)
          ] *
            dt
        )
      
      
      # --------------------------------------------------------
      # Boundary for R at infection age zero:
      #
      # R_kj(t,0)
      #
      # =
      #
      # beta (j-1) PSS_kj Lambda_j / Lambda_k
      #
      # This is well defined in the nonlinear calculation.
      # --------------------------------------------------------
      
      R_new[
        ik,
        jj,
        1
      ] <-
        beta *
        (j - 1) *
        PSS_new[
          ik,
          jj
        ] *
        Lambda_new[jj] /
        max(
          Lambda_new[ik],
          1e-14
        )
      
    }

    
  }
  
  
  # ============================================================
  # J. UPDATE P_{S_j}^{I_k}
  # ============================================================
  
  PIS_new <- array(
    0,
    dim = c(K, K, Na)
  )
  
  
  for (ik in 1:K) {
    
    k <- degrees[ik]
    
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      
      rate_PIS <-
        -beta -
        beta *
        (j - 1) *
        Lambda[jj] +
        delta_p *
        PT_total[
          ik,
          
        ]
      
      
      PIS_new[
        ik,
        jj,
        2:Na
      ] <-
        PIS[
          ik,
          jj,
          1:(Na - 1)
        ] *
        exp(
          rate_PIS[
            1:(Na - 1)
          ] *
            dt
        )
      
      
      # Infection-age boundary
      
      PIS_new[
        ik,
        jj,
        1
      ] <-
        ((k - 1) / k) *
        PSS_new[
          ik,
          jj
        ]
      
    }
    
  }
  
  
  PIS_new <- pmax(
    PIS_new,
    0
  )
  
  
  # ============================================================
  # K. UPDATE P_{I_j}^{I_k}(t,a1,a2)
  # ============================================================
  
  PII_new <- array(
    0,
    dim = c(
      K,
      K,
      Na,
      Na
    )
  )
  
  
  for (ik in 1:K) {
    
    k <- degrees[ik]
    
    
    for (jj in 1:K) {
      
      j <- degrees[jj]
      
      
      # --------------------------------------------------------
      # Interior two-age characteristic equation
      # --------------------------------------------------------
      
      rate_matrix <-
        outer(
          PT_total[
            ik,
            
          ],
          PT_total[
            jj,
            
          ],
          FUN = function(
    PT_k,
    PT_j
          ) {
            
            -(gamma + alpha) -
              delta_p *
              (j - 1) *
              PT_j +
              delta_p *
              PT_k
            
          }
        )
      
      
      PII_new[
        ik,
        jj,
        2:Na,
        2:Na
      ] <-
        PII[
          ik,
          jj,
          1:(Na - 1),
          1:(Na - 1)
        ] *
        exp(
          rate_matrix[
            1:(Na - 1),
            1:(Na - 1)
          ] *
            dt
        )
      
      
      # --------------------------------------------------------
      # BOUNDARY 1 -- GOLD-STANDARD NORMALIZATION
      #
      # P_{I_j}^{I_k}(t,0,a2)
      #
      #  =
      #
      # [ (k-1) Lambda_k + 1/k ]
      # *
      # R_kj(t,a2)
      #
      # obtained from:
      # P_{I_j}^{S_k}
      # * (k - 1 + 1/(k Lambda_k))
     
      # --------------------------------------------------------
      
      PII_new[
        ik,
        jj,
        1,
        
      ] <-
        (
          (1 + (k - 1) * Lambda_new[ik]) / k
        ) *
        R_new[
          ik,
          jj,

        ]
      
      
      # --------------------------------------------------------
      # BOUNDARY 2
      #
      # P_{I_j}^{I_k}(t,a1,0)
      #
      # =
      #
      # beta [1+(j-1)Lambda_j]
      # P_{S_j}^{I_k}(t,a1)
      # --------------------------------------------------------
      
      PII_new[
        ik,
        jj,
        ,
        1
      ] <-
        beta *
        (
          1 +
            (j - 1) *
            Lambda_new[jj]
        ) *
        PIS_new[
          ik,
          jj,
          
        ]
      
    }
    
  }
  
  
  PII_new <- pmax(
    PII_new,
    0
  )
  
  
  # ============================================================
  # L. UPDATE P_{T_j}^{I_k}: GOLD-STANDARD EQUATION
  #
  # (d_t+d_a) P_{T_j}^{I_k}
  # = -delta P_{T_j}^{I_k}
  #   + alpha * integral P_{I_j}^{I_k}(a1,a2) da2
  #   + delta*p*(j-1) * integral
  #       P_{I_j}^{I_k}(a1,a2) P_T^{I_j}(a2) da2
  #   + delta*p P_{T_j}^{I_k} P_T^{I_k}.
  #
  # No (2k-1) closure and no probability projection.
  # ============================================================

  PTI_new <- array(0, dim = c(K, K, Na))

  for (ik in 1:K) {
    k <- degrees[ik]

    for (jj in 1:K) {
      j <- degrees[jj]

      W_old <- W[ik, jj, ]

      tracing_creation <- numeric(Na)
      for (a1 in 1:Na) {
        tracing_creation[a1] <-
          simpson_age(
            PII[ik, jj, a1, ] *
              PT_total[jj, ]
          )
      }

      # dP/dt = coeff*P + source along each characteristic.
      coeff <-
        -delta +
        delta_p * PT_total[ik, ]

      source <-
        alpha * W_old +
        delta_p * (j - 1) * tracing_creation

      c_old <- coeff[1:(Na - 1)]
      b_old <- source[1:(Na - 1)]

      phi <- ifelse(
        abs(c_old) > 1e-12,
        (exp(c_old * dt) - 1) / c_old,
        dt
      )

      PTI_new[ik, jj, 2:Na] <-
        exp(c_old * dt) *
          PTI[ik, jj, 1:(Na - 1)] +
        phi * b_old

      # P_{T_j}^{I_k}(t,0) = (k-1)/k P_{T_j}^{S_k}(t)
      PTI_new[ik, jj, 1] <-
        ((k - 1) / k) * PST_new[ik, jj]
    }
  }

  # Numerical roundoff guard only; no rescaling/projection.
  PTI_new <- pmax(PTI_new, 0)

  # ============================================================
  # N. NUMERICAL SAFETY TEST
  # ============================================================
  
  finite_test <-
    all(
      is.finite(
        I_new
      )
    ) &&
    all(
      is.finite(
        Lambda_new
      )
    ) &&
    all(
      is.finite(
        R_new
      )
    ) &&
    all(
      is.finite(
        PIS_new
      )
    ) &&
    all(
      is.finite(
        PII_new
      )
    ) &&
    all(
      is.finite(
        PTI_new
      )
    )
  
  
  if (!finite_test) {
    
    cat("\n")
    cat("NON-FINITE VALUE DETECTED\n")
    cat(
      "time =",
      times[n],
      "\n"
    )
    
    stop(
      "Numerical instability detected."
    )
    
  }
  
  
  # ============================================================
  # O. UPDATE SYSTEM STATE
  # ============================================================
  
  S <- S_new
  
  I_age <- I_new
  
  PSS <- PSS_new
  
  PST <- PST_new
  
  Lambda <- Lambda_new
  
  R <- R_new
  
  PIS <- PIS_new
  
  PII <- PII_new
  
  PTI <- PTI_new
  
  
  Lambda_history[
    n + 1,
    
  ] <- Lambda
  
  
  # ============================================================
  # P. TOTAL INFECTIOUS POPULATION
  # ============================================================
  
  I_zero <-
    I_zero_0 *
    exp(
      -(alpha + gamma) *
        times[n + 1]
    )
  
  
  I_total[
    n + 1
  ] <-
    sum(
      apply(
        I_age,
        1,
        simpson_age
      )
    ) +
    I_zero
  
}


# ==============================================================
# 14. FINISHED
# ==============================================================

cat("\n")
cat("============================================\n")
cat("PDE simulation completed successfully!\n")
cat("============================================\n\n")


# ==============================================================
# 15. SAVE PDE RESULTS
# ==============================================================

pde_data <- data.frame(
  t = times,
  I_PDE = I_total
)


write.csv(
  pde_data,
  file = "PDE_It_gold_standard_fast.csv",
  row.names = FALSE
)


# ==============================================================
# 16. PDE DIAGNOSTICS
# ==============================================================

cat(
  "Initial I(t) =",
  I_total[1],
  "\n"
)

cat(
  "Peak PDE I(t) =",
  max(
    I_total,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "PDE peak time =",
  times[
    which.max(
      I_total
    )
  ],
  "\n"
)

cat(
  "Final PDE I(t) =",
  tail(
    I_total,
    1
  ),
  "\n"
)


# ==============================================================
# 17. PDE-ONLY PLOT
# ==============================================================

plot(
  times,
  I_total,
  type = "l",
  lwd = 2,
  xlab = "Time",
  ylab = "Number infectious",
  log = "y",
  main = 
    "Nonlinear Age-Structured PDE Contact-Tracing Model"
)


# ==============================================================
# 18. COMPARE AGAINST PROFESSOR'S ABM
# ==============================================================

if (
  file.exists(
    "It_p030.csv"
  )
) {
  
  
  abm <-
    read.csv(
      "It_p030.csv"
    )
  
  
  I_ABM <-
    rowMeans(
      abm[
        ,
        -1,
        drop = FALSE
      ]
    )
  
  
  cat("\n")
  cat(
    "Peak ABM I(t) =",
    max(
      I_ABM
    ),
    "\n"
  )
  
  cat(
    "ABM peak time =",
    abm$t[
      which.max(
        I_ABM
      )
    ],
    "\n"
  )
  
  
  ymax <-
    max(
      c(
        I_ABM,
        I_total
      ),
      na.rm = TRUE
    )
  
  
  plot(
    abm$t,
    I_ABM,
    type = "l",
    lwd = 2,
    xlab = "Time",
    ylab = "Number infectious",
    ylim = c(
      0,
      1.05 * ymax
    ),
    main =
      "PDE versus Agent-Based Simulation"
  )
  
  
  lines(
    times,
    I_total,
    lwd = 2,
    lty = 2
  )
  
  
  legend(
    "topright",
    legend = c(
      "Agent-based model",
      "PDE model"
    ),
    lty = c(
      1,
      2
    ),
    lwd = c(
      2,
      2
    ),
    bty = "n"
  )
  
  
} else {
  
  
  cat("\n")
  cat(
    "It.csv was not found.\n"
  )
  
  cat(
    "The PDE itself completed successfully.\n"
  )
  
}


# ==============================================================
# 19. FINAL MESSAGE
# ==============================================================

cat("\n")
cat(
  "Saved PDE results to PDE_It_gold_standard_fast.csv\n"
)

cat(
  "DONE.\n"
)