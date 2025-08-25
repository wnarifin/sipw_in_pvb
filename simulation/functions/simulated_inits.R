# Inits ====

# General Init ====
seednum = 3209673  # random digit from electricity bill
# combine harel 2006, rochani 2015 settings
# n = c(50, 100, 200, 500, 1000)  # to test 50, 100, 200, 500 and 1000
# n1 = c(100, 200, 500, 1000)  # to test 100, 200, 500 and 1000, larger size for with covar, mnar
# n2 = c(200, 500, 1000)  # to test 200, 500 and 1000, larger size for with covar, mnar
n = c(200, 500, 1000)  # to test 200, 500 and 1000, larger size for with covar, mnar
# B = 2000  # number of different samples/simulated replicates -- rochani2015
# m = 500  # number of imputations
B = 500  # number of different samples/simulated replicates -- unal&burgut2014
# B = 1000  # number of different samples/simulated replicates -- he & mcdermott 2012
m = 100  # number of imputations
# have to settle for a lower number, else will take long time for imputations + bootstrap
b = 1000  # number of bootstraps -- woodward2014, unal & burgutt 2014 (=1000)
# for CI estimates, use b=999
# b = 100  # he & mcdermott 2012

# Init pop values ====
# possible values
# sn = c(0.2, 0.3, 0.6, 0.75, 0.9, 0.95)  # add 0.2, He & McDermott
# sp = c(0.6, 0.75, 0.8, 0.9, 0.95)
# p = c(0.03, 0.1, 0.4, 0.5)  # add 0.1, He & McDermott
# must consider reasonable combinations that are distinct
# not too many combination, else not manageable

# # Init 1a, Harel & Zhou 1, default
# p_pop = 0.4; Sn_pop = 0.9; Sp_pop = 0.8
# 
# # Init 1b, Harel & Zhou 1
# p_pop = 0.4; Sn_pop = 0.95; Sp_pop = 0.8
# 
# # Init 2, Harel & Zhou 2
# p_pop = 0.03; Sn_pop = 0.3; Sp_pop = 0.9

# High Sn Sp
Sn_pop_list = c(0.3, 0.6, 0.9)  # incr by 0.3
Sp_pop_list = c(0.6, 0.9)  # incr by 0.3
p_pop_list = c(0.05, 0.1, 0.4)  # reasonable values

# For each p_pop_list, combinations of Sn_pop with Sp_pop
SnSp_pop_df = data.frame(Sn = rep(Sn_pop_list, each = length(Sp_pop_list)), 
                         Sp = rep(Sp_pop_list, length(Sn_pop_list)))
SnSp_pop_df

# Init for MAR ====
# Pr verification when T=t

# # Init 1, Harel & Zhou 1, default
# p_sel_t1 = 0.8  # P(V=1|T=1)
# p_sel_t0 = 0.4  # P(V=1|T=0)
# 
# # Init 2, Harel & Zhou 2
# p_sel_t1 = 0.55  # P(V=1|T=1)
# # p_sel_t0 = 0.06  # P(V=1|T=0)
# p_sel_t0 = 0.05  # P(V=1|T=0)  # modified

# DF for p_sel_t
p_sel_t1_list = c(0.55, 0.8)
p_sel_t0_list = c(0.05, 0.4)
p_sel_t_df = data.frame(t1 = p_sel_t1_list, t0 = p_sel_t0_list)
p_sel_t_df

# # Set Ps for simple 2x2 table ====
# 
# # multinomial distribution
# # cell Ps calculated based on Bayesian rules
# # p_td
# # p_+1 = p, p_+0 = 1-p
# # p_1+ = t+, p_0+ = t-
# p11 = Sn_pop * p_pop  # s1: P(T=1, D=1) = P(T=1|D=1)P(D=1)
# p01 = (1-Sn_pop) * p_pop # s0: P(T=0, D=1) = P(T=0|D=1)P(D=1)
# p10 = (1-Sp_pop) * (1-p_pop)  # r1: P(T=1, D=0) = P(T=1|D=0)P(D=0)
# p00 = Sp_pop * (1-p_pop)  # r0: P(T=0, D=0) = P(T=0|D=0)P(D=0)
# 
# # Define data generating distribution ====
# # cells by col, i.e. by D=d
# def = defData(varname = "td", formula = "p11;p01;p10;p00", dist = "categorical")

