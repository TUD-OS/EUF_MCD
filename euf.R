#Calculate idle model ...
#Use visreg (package visreg) to visualize regression
required.packages <- c("crayon","RSQLite","proto","gsubfn","readr","sqldf")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for (p in required.packages) {
  library(p,character.only =TRUE)
}

print_eqn <- function(sum, prefix = "") {
  for (r in rownames(coef(sum))[2:length(rownames(coef(sum)))]) {
    cat(coef(sum)[r,"Estimate"])
    for (prod in strsplit(r,":")[[1]]) {
      cat(" * ")
      #cat("\nDEBUG: ",r,"=>",prod,"\n")
      if (substr(prod,1,4) == "poly") {
        elem <- strsplit(sub("poly\\(([^,]+), [^\\)]+\\)([0-9]*)","\\1,\\2",prod, perl=TRUE),",")[[1]]
        cat(prefix,elem[1],sep="")
        cat(ifelse(elem[2] > 1,paste("**",elem[2],sep=""),""))
      } else {
        cat(prefix,prod,sep="")
      }
    }
    cat(" + ",sep="")
  }
  cat(coef(sum)["(Intercept)","Estimate"])
  str = paste(" [R² =",(sum$adj.r.squared),"]\n")
  if (sum$adj.r.squared > 0.9) {
    cat(green(str))
  } else if (sum$adj.r.squared > 0.8) {
    cat(yellow(str))
  } else if (sum$adj.r.squared > 0.7) {
    cat(str)
  } else if (sum$adj.r.squared > 0.6) {
    cat(magenta(str))
  } else {
    cat(red(str))
  }
}

#solve_eqn <- function(sum,values) {
solve_eqn <- function(sum,...) {
  tot_sum = 0
  for (r in rownames(coef(sum))[2:length(rownames(coef(sum)))]) {
    curval <- coef(sum)[r,"Estimate"]
    for (prod in strsplit(r,":")[[1]]) {
      if (substr(prod,1,4) == "poly") {
        elem <- strsplit(sub("poly\\(([^,]+), [^\\)]+\\)([0-9]*)","\\1,\\2",prod, perl=TRUE),",")[[1]]
        #curval <- curval * values[[elem[1]]] ^ as.numeric(elem[2])
        curval <- curval * (list(...)[[elem[1]]] ^ as.numeric(elem[2]))
      } else {
        #curval <- curval * values[[elem[1]]]
        curval <- curval * list(...)[[prod]]
      }
    }
    tot_sum <- tot_sum + curval
  }
  tot_sum <- tot_sum + coef(sum)["(Intercept)","Estimate"]
  return(tot_sum)
}

idle_new <- read_delim("idle_new.csv", ";", escape_double = FALSE, trim_ws = TRUE,col_types = cols())
idle_new = within(idle_new, {
  Frequency <- Frequency/1000; #Freq in MHz
  #Divide by ten because there were ten rounds! E is measured over all
  P_arm = E_arm/time/10
  P_mem = E_mem/time/10
  P_g3d = E_g3d/time/10
  P_kfc = E_kfc/time/10
  P_ext = E_ext/time/10
  kIPS = instructions_k/time
  uIPS = instructions_u/time
  P_Diff = P_ext - P_arm - P_g3d - P_kfc - P_mem;
})

m_ARM <- lm(P_arm~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "big_global_idle_bigNet"')) 
m_KFC <- lm(P_kfc~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "little_global_idle"'))
m_kIPS_Idle_Big <- lm(kIPS~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "big_global_idle"'))
m_uIPS_Idle_Big <- lm(uIPS~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "big_global_idle"'))
m_kIPS_Idle_Little <- lm(kIPS~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "little_global_idle_bigNet"'))
m_uIPS_Idle_Little <- lm(uIPS~poly(Frequency,2,raw=TRUE),data=sqldf('select * from idle_new where Config = "little_global_idle_bigNet"'))
                                        
sm_ARM <- summary(m_ARM)
sm_KFC <- summary(m_KFC)
sm_kIPS_Idle_Big <- summary(m_kIPS_Idle_Big)
sm_uIPS_Idle_Big <- summary(m_uIPS_Idle_Big)
sm_kIPS_Idle_Little <- summary(m_kIPS_Idle_Little)
sm_uIPS_Idle_Little <- summary(m_uIPS_Idle_Little)

P_Idle_ARM <- function(freq){
  return (sm_ARM$coefficients[2]*freq+sm_ARM$coefficients[3]*freq*freq+sm_ARM$coefficients[1])
}

P_Idle_KFC <- function(freq){
  return (sm_KFC$coefficients[2]*freq+sm_KFC$coefficients[3]*freq*freq+sm_KFC$coefficients[1])
}

kInsns <- function(freq,conf){
  if (startsWith(conf,"big")[1]) {
    return(sqldf(paste('select * from idle_new where Frequency == "',freq,'" and Config = "',conf,'_global_idle"',sep=""))$kIPS)
  } else {
    return(sqldf(paste('select * from idle_new where Frequency == "',freq,'" and Config = "',conf,'_global_idle_bigNet"', sep = ""))$kIPS)
  }
}

uInsns <- function(freq,conf){
  if (startsWith(conf,"big")[1]) {
    return(sqldf(paste('select * from idle_new where Frequency == "',freq,'" and Config = "',conf,'_global_idle"',sep = ""))$uIPS)
  } else {
    return(sqldf(paste('select * from idle_new where Frequency == "',freq,'" and Config = "',conf,'_global_idle_bigNet"',sep=""))$uIPS)
  }
}

getMeasured <- function(Frequency,CacheSize,CacheThreads,field) {
    return (sqldf(paste('select ',field,' from bench_orig where Frequency=',Frequency,' and CacheThreads=',CacheThreads,' and CacheSize=',CacheSize,' and Config="big"',sep=""))[[field]])
}

#cat("Hello: ",uInsns(200000,"big"),"\n")

P_Idle_MEM <- mean(idle_new$P_mem)
P_Idle_G3D <- mean(idle_new$P_g3d)
P_Idle_Rest <- mean(idle_new$P_Diff)

cat("P_Idle_ARM = "); print_eqn(sm_ARM,prefix="Freq_");
cat("P_Idle_KFC = "); print_eqn(sm_KFC,prefix="Freq_");
cat("P_Idle_MEM =",P_Idle_MEM,"W\n")
cat("P_Idle_G3D =",P_Idle_G3D,"W\n")
cat("P_Idle_Rest =",P_Idle_Rest,"W\n")
#cat("kIPC_Idle_Big = "); print_eqn(sm_kIPS_Idle_Big)
#cat("uIPC_Idle_Big = "); print_eqn(sm_uIPS_Idle_Big)
#cat("kIPC_Idle_Little = "); print_eqn(sm_kIPS_Idle_Little)
#cat("uIPC_Idle_Little = "); print_eqn(sm_uIPS_Idle_Little)

#TODO: Replace bak by pp
bench <- read_delim("bench.csv", ";", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
bench_orig <- within(bench, {
  requests = Gets + Sets
  E_total = E_arm + E_g3d + E_mem + E_kfc
  P_ext = E_ext/t_wall
  rps = requests / t_wall
  kIPR = instructions_k / requests
  uIPR = instructions_u / requests
  IPR = uIPR + kIPR
  instructions = instructions_k + instructions_u
  memoryHeaviness = (`LLC-load-misses` + `LLC-store-misses`)/instructions
})
bench <- sqldf("select * from bench where CacheSize <> 128")
bench = within(bench, {
  #Frequency = Frequency*1000
  requests = Gets + Sets
  E_total = E_arm + E_g3d + E_mem + E_kfc
  P = E_total / t_wall
  P_ext = E_ext/t_wall
  instructions_k_adj = instructions_k - kInsns(Frequency,Config)*t_wall
  instructions_u_adj = instructions_u - uInsns(Frequency,Config)*t_wall
  kIPR = instructions_k / requests
  uIPR = instructions_u / requests
  kIPC = instructions_k / t_wall / (Frequency * 1E6)
  uIPC = instructions_u / t_wall / (Frequency * 1E6)
  EperRequest = E_total / requests
  rps = requests / t_wall
  P_mem = E_mem / t_wall
  P_arm = E_arm / t_wall
  P_g3d = E_g3d / t_wall
  P_kfc = E_kfc / t_wall
  P_dyn_arm = P_arm - ifelse(Config == 'big', P_Idle_ARM(Frequency), P_arm)
  P_dyn_kfc = P_kfc - ifelse(Config == 'little', P_Idle_KFC(Frequency), P_kfc)
  P_dyn_mem = P_mem - P_Idle_MEM
  P_dyn_g3d = P_g3d - P_Idle_G3D
  P_Rest = P_ext - P_arm - P_g3d - P_mem - P_kfc
  P_dyn_Rest = P_Rest - P_Idle_Rest
  IPC = uIPC + kIPC
  IPR = uIPR + kIPR
  instructions = instructions_k + instructions_u
  instructions_adj = instructions_k_adj + instructions_u_adj
  memoryHeaviness = (`LLC-load-misses` + `LLC-store-misses`)/instructions
})
m_simple_P <- lm(P_ext ~ poly(Frequency,3,raw=TRUE)+poly(CacheThreads,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
m_simple_rps <- lm(rps ~ poly(Frequency,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))

#bench <- sqldf('select * from bench where instructions_u > 5e9')
#m_MEM_Dyn <- lm(P_dyn_mem~poly(memoryHeaviness,2,raw=TRUE)+poly(Frequency,2,raw=TRUE),data=bench)
m_MEM_Dyn <- lm(P_dyn_mem~poly(memoryHeaviness,2,raw=TRUE)+poly(Frequency,2,raw=TRUE)+poly(IPC,2,raw=TRUE),data=bench)
m_MEM_Dyn_Big <- lm(P_dyn_mem~poly(Frequency,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
m_MEM_Dyn_Little <- lm(P_dyn_mem~poly(Frequency,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
#m_ARM_Dyn <- lm(P_dyn_arm~poly(IPC,2,raw=TRUE)+poly(Frequency,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
m_ARM_Dyn <- lm(P_dyn_arm~poly(Frequency,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE)+poly(IPC,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
#m_KFC_Dyn <- lm(P_dyn_kfc~poly(IPC,2,raw=TRUE)+poly(Frequency,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
m_KFC_Dyn <- lm(P_dyn_kfc~poly(Frequency,2,raw=TRUE)+CacheThreads,data=sqldf('select * from bench where Config == "little"'))
m_Rest_Dyn <- lm(P_dyn_Rest~poly(Frequency,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))

sm_simple_P <- summary(m_simple_P)
sm_simple_rps <- summary(m_simple_rps)
sm_MEM_Dyn <- summary(m_MEM_Dyn)
sm_MEM_Dyn_Big <- summary(m_MEM_Dyn_Big)
sm_MEM_Dyn_Little <- summary(m_MEM_Dyn_Little)
sm_ARM_Dyn <- summary(m_ARM_Dyn)
sm_KFC_Dyn <- summary(m_KFC_Dyn)
sm_Rest_Dyn <- summary(m_Rest_Dyn)

cat("P_Dyn_ARM = "); print_eqn(sm_ARM_Dyn)
cat("P_Dyn_KFC = "); print_eqn(sm_KFC_Dyn)
cat("P_Dyn_MEM = "); print_eqn(sm_MEM_Dyn)
#cat("P_Dyn_MEM_Big = "); print_eqn(sm_MEM_Dyn_Big)
#cat("P_Dyn_MEM_Little = "); print_eqn(sm_MEM_Dyn_Little)
cat("P_Dyn_Rest = "); print_eqn(sm_Rest_Dyn)
cat("P_simple = "); print_eqn(sm_simple_P)
cat("rps_simple = "); print_eqn(sm_simple_rps)


#Application Model
#m_MCD_IPC <- lm(IPC~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=bench)
#m_MCD_IPC_Big <- lm(IPC~poly(rps,2,raw=TRUE)*poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
#m_MCD_IPC_Big <- lm(IPC~poly(rps,2,raw=TRUE)*IPR,data=sqldf('select * from bench where Config == "big"'))
m_MCD_IPC_Big2 <- lm(IPC~poly(Frequency,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE)+memoryHeaviness,data=sqldf('select * from bench where Config == "big"'))
#m_MCD_IPC_Little <- lm(IPC~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
#m_MCD_IPC_Little <- lm(IPC~poly(rps,2,raw=TRUE)*IPR,data=sqldf('select * from bench where Config == "little"'))

#m_MCD_kIPC <- lm(kIPC~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=bench)
#m_MCD_kIPC_Big <- lm(kIPC~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
#m_MCD_kIPC_Little <- lm(kIPC~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))


#sm_MCD_IPC <- summary(m_MCD_IPC)
#sm_MCD_IPC_Big <- summary(m_MCD_IPC_Big)
sm_MCD_IPC_Big2 <- summary(m_MCD_IPC_Big2)
#sm_MCD_IPC_Little <- summary(m_MCD_IPC_Little)
#m_MCD_IPR <- lm(IPR~poly(rps,2,raw=TRUE)+poly(CacheSize,2,raw=TRUE),data=bench)
m_MCD_IPR_Big <- lm(IPR~poly(CacheSize,2,raw=TRUE)+CacheThreads,data=sqldf('select * from bench where Config == "big"'))
m_MCD_IPR_Big_orig <- lm(IPR~poly(CacheSize,2,raw=TRUE)+CacheThreads,data=sqldf('select * from bench_orig where Config == "big"'))
#m_MCD_IPR_Little <- lm(IPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
#m_MCD_kIPR <- lm(kIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=bench)
#m_MCD_kIPR_Big <- lm(kIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
#m_MCD_kIPR_Little <- lm(kIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
#m_MCD_uIPR <- lm(uIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=bench)
#m_MCD_uIPR_Big <- lm(uIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
#m_MCD_uIPR_Little <- lm(uIPR~poly(CacheSize,2,raw=TRUE)+poly(CacheThreads,2,raw=TRUE),data=sqldf('select * from bench where Config == "little"'))
m_MCD_memoryHeaviness_Big <- lm(memoryHeaviness~poly(CacheSize,2,raw=TRUE)+poly(IPR,2,raw=TRUE),data=sqldf('select * from bench where Config == "big"'))
m_MCD_memoryHeaviness_Big_orig <- lm(memoryHeaviness~poly(CacheSize,2,raw=TRUE)+poly(IPR,2,raw=TRUE),data=sqldf('select * from bench_orig where Config == "big"'))
#m_MCD_memoryHeaviness_Big <- lm(memoryHeaviness~IPR,data=sqldf('select * from bench where Config == "big"'))
#sm_MCD_IPR <- summary(m_MCD_IPR)
sm_MCD_IPR_Big <- summary(m_MCD_IPR_Big)
sm_MCD_IPR_Big_orig <- summary(m_MCD_IPR_Big_orig)
#sm_MCD_IPR_Little <- summary(m_MCD_IPR_Little)
#sm_MCD_kIPR <- summary(m_MCD_kIPR)
#sm_MCD_kIPR_Big <- summary(m_MCD_kIPR_Big)
#sm_MCD_kIPR_Little <- summary(m_MCD_kIPR_Little)
#sm_MCD_uIPR <- summary(m_MCD_uIPR)
#sm_MCD_uIPR_Big <- summary(m_MCD_uIPR_Big)
#sm_MCD_uIPR_Little <- summary(m_MCD_uIPR_Little)
sm_MCD_memoryHeaviness_Big <- summary(m_MCD_memoryHeaviness_Big)
sm_MCD_memoryHeaviness_Big_orig <- summary(m_MCD_memoryHeaviness_Big_orig)

#cat("IPC_MCD = "); print_eqn(sm_MCD_IPC)
#cat("IPC_MCD_Big = "); print_eqn(sm_MCD_IPC_Big)
cat("IPC_MCD_Big2 = "); print_eqn(sm_MCD_IPC_Big2)
#cat("IPC_MCD_Little = "); print_eqn(sm_MCD_IPC_Little)
#cat("IPR_MCD = "); print_eqn(sm_MCD_IPR)
cat("IPR_MCD_Big = "); print_eqn(sm_MCD_IPR_Big)
#cat("IPR_MCD_Little = "); print_eqn(sm_MCD_IPR_Little)
#cat("kIPR_MCD = "); print_eqn(sm_MCD_kIPR)
#cat("kIPR_MCD_Big = "); print_eqn(sm_MCD_kIPR_Big)
#cat("kIPR_MCD_Little = "); print_eqn(sm_MCD_kIPR_Little)
#cat("uIPR_MCD = "); print_eqn(sm_MCD_uIPR)
#cat("uIPR_MCD_Big = "); print_eqn(sm_MCD_uIPR_Big)
#cat("uIPR_MCD_Little = "); print_eqn(sm_MCD_uIPR_Little)
cat("memoryHeaviness_Big = "); print_eqn(sm_MCD_memoryHeaviness_Big)
write_delim(bench,"bench.result.csv",";")

model_light = expand.grid(Frequency = seq(200,2000,100), Cores = c(1,2), CacheSize = c(32,64,256))
model_light = within(model_light, {
    IPR = solve_eqn(sm_MCD_IPR_Big,CacheSize = CacheSize, CacheThreads = Cores,Frequency=Frequency)
    IPR_orig = solve_eqn(sm_MCD_IPR_Big_orig,CacheSize = CacheSize, CacheThreads = Cores,Frequency=Frequency)
    memory_heaviness = solve_eqn(sm_MCD_memoryHeaviness_Big, CacheSize = CacheSize, CacheThreads = Cores,IPR=IPR)
    memory_heaviness_orig = solve_eqn(sm_MCD_memoryHeaviness_Big_orig, CacheSize = CacheSize, CacheThreads = Cores,IPR=IPR_orig)
    IPC = solve_eqn(sm_MCD_IPC_Big2,Frequency = Frequency, CacheThreads = Cores,CacheSize=CacheSize, memoryHeaviness = memory_heaviness)
    IPC_orig = solve_eqn(sm_MCD_IPC_Big2,Frequency = Frequency, CacheThreads = Cores,CacheSize=CacheSize, memoryHeaviness = memory_heaviness_orig)
    rps = Frequency * 1000000 / IPR * IPC
    rps_orig = Frequency * 1000000 / IPR_orig * IPC_orig
    P_idle_arm = solve_eqn(sm_ARM, Frequency = Frequency)
    P_idle_kfc = solve_eqn(sm_KFC, Frequency = Frequency)
    P_idle = P_idle_arm + P_idle_kfc + P_Idle_MEM + P_Idle_G3D + P_Idle_Rest
    P_dyn_arm = solve_eqn(sm_ARM_Dyn,Frequency = Frequency, CacheThreads = Cores,IPC=IPC)
    P_dyn_arm_orig = solve_eqn(sm_ARM_Dyn,Frequency = Frequency, CacheThreads = Cores,IPC=IPC_orig)
    P_dyn_mem = solve_eqn(sm_MEM_Dyn,Frequency = Frequency, IPC = IPC,memoryHeaviness = memory_heaviness, IPR=IPR)
    P_dyn_mem_orig = solve_eqn(sm_MEM_Dyn,Frequency = Frequency, IPC = IPC_orig,memoryHeaviness = memory_heaviness_orig, IPR=IPR_orig)
    P_dyn_rest = solve_eqn(sm_Rest_Dyn,Frequency = Frequency)
    P_dyn = P_dyn_arm + P_dyn_mem + P_dyn_rest
    P_dyn_orig = P_dyn_arm_orig + P_dyn_mem_orig + P_dyn_rest
    P = P_dyn + P_idle
    P_orig = P_dyn_orig + P_idle
})

model = expand.grid(Frequency = seq(200,2000,100), Cores = c(1,2), CacheSize = c(32,64,128,256))
model = within(model, {
    IPR = solve_eqn(sm_MCD_IPR_Big,CacheSize = CacheSize, CacheThreads = Cores,Frequency=Frequency)
    IPR_orig = solve_eqn(sm_MCD_IPR_Big_orig,CacheSize = CacheSize, CacheThreads = Cores,Frequency=Frequency)
    memory_heaviness = solve_eqn(sm_MCD_memoryHeaviness_Big, CacheSize = CacheSize, CacheThreads = Cores,IPR=IPR)
    memory_heaviness_orig = solve_eqn(sm_MCD_memoryHeaviness_Big_orig, CacheSize = CacheSize, CacheThreads = Cores,IPR=IPR_orig)
    IPC = solve_eqn(sm_MCD_IPC_Big2,Frequency = Frequency, CacheThreads = Cores,CacheSize=CacheSize, memoryHeaviness = memory_heaviness)
    IPC_orig = solve_eqn(sm_MCD_IPC_Big2,Frequency = Frequency, CacheThreads = Cores,CacheSize=CacheSize, memoryHeaviness = memory_heaviness_orig)
    rps = Frequency * 1000000 / IPR * IPC
    rps_orig = Frequency * 1000000 / IPR_orig * IPC_orig
    P_idle_arm = solve_eqn(sm_ARM, Frequency = Frequency)
    P_idle_kfc = solve_eqn(sm_KFC, Frequency = Frequency)
    P_idle = P_idle_arm + P_idle_kfc + P_Idle_MEM + P_Idle_G3D + P_Idle_Rest
    P_dyn_arm = solve_eqn(sm_ARM_Dyn,Frequency = Frequency, CacheThreads = Cores,IPC=IPC)
    P_dyn_arm_orig = solve_eqn(sm_ARM_Dyn,Frequency = Frequency, CacheThreads = Cores,IPC=IPC_orig)
    P_dyn_mem = solve_eqn(sm_MEM_Dyn,Frequency = Frequency, IPC = IPC,memoryHeaviness = memory_heaviness, IPR=IPR)
    P_dyn_mem_orig = solve_eqn(sm_MEM_Dyn,Frequency = Frequency, IPC = IPC_orig,memoryHeaviness = memory_heaviness_orig, IPR=IPR_orig)
    P_dyn_rest = solve_eqn(sm_Rest_Dyn,Frequency = Frequency)
    P_dyn = P_dyn_arm + P_dyn_mem + P_dyn_rest
    P_dyn_orig = P_dyn_arm_orig + P_dyn_mem_orig + P_dyn_rest
    P = P_dyn + P_idle
    P_orig = P_dyn_orig + P_idle
})

model_mono_light = expand.grid(Frequency = seq(200,2000,100), Cores = c(1,2), CacheSize = c(32,64,256))
model_mono_light = within(model_mono_light, {
    rps = solve_eqn(sm_simple_rps,CacheSize = CacheSize, Frequency = Frequency, CacheThreads = Cores)
    P = solve_eqn(sm_simple_P,CacheSize = CacheSize, Frequency = Frequency, CacheThreads = Cores)
})
model_mono = expand.grid(Frequency = seq(200,2000,100), Cores = c(1,2), CacheSize = c(32,64,128,256))
model_mono = within(model_mono, {
    rps = solve_eqn(sm_simple_rps,CacheSize = CacheSize, Frequency = Frequency, CacheThreads = Cores)
    P = solve_eqn(sm_simple_P,CacheSize = CacheSize, Frequency = Frequency, CacheThreads = Cores)
})

model$measured_P <- apply(model,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'P_ext'))
model_light$measured_P <- apply(model_light,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'P_ext'))
model_mono$measured_P <- apply(model_mono,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'P_ext'))
model_mono_light$measured_P <- apply(model_mono_light,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'P_ext'))
model$measured_rps <- apply(model,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'rps'))
model_light$measured_rps <- apply(model_light,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'rps'))
model_mono$measured_rps <- apply(model_mono,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'rps'))
model_mono_light$measured_rps <- apply(model_mono_light,1,function(x) getMeasured(Frequency=x[['Frequency']], CacheThreads=x[['Cores']], CacheSize=x[['CacheSize']], 'rps'))

model = within(model, {
    delta_P = measured_P - P
    delta_rps = measured_rps - rps
    delta_P_rel = delta_P / measured_P
    delta_rps_rel = delta_rps / measured_rps
    #Recalibrated
    delta_P_orig = measured_P - P_orig
    delta_rps_orig = measured_rps - rps_orig
    delta_P_rel_orig = delta_P_orig / measured_P
    delta_rps_rel_orig = delta_rps_orig / rps_orig
})
model_light = within(model_light, {
    delta_P = measured_P - P
    delta_rps = measured_rps - rps
    delta_P_rel = delta_P / P
    delta_rps_rel = delta_rps / rps
})
model_mono_light = within(model_mono_light, {
    delta_P = measured_P - P
    delta_rps = measured_rps - rps
    delta_P_rel = delta_P / P
    delta_rps_rel = delta_rps / rps
})
model_mono = within(model_mono, {
    delta_P = measured_P - P
    delta_rps = measured_rps - rps
    delta_P_rel = delta_P / P
    delta_rps_rel = delta_rps / rps
})



write_delim(model,"bench.modeled.csv",";")
write_delim(model_light,"bench.modeled_light.csv",";")
write_delim(model_mono,"bench.modeled_mono.csv",";")
write_delim(model_mono_light,"bench.modeled_mono_light.csv",";")

