# wf function ----------------------------------------------

wf <- function(npop=5,ngen=500,nall=100,init=1,mu=0)
  # this is a comment. Comments inform users; they are ignored by R
{
  freq = matrix(init,npop,ngen) # create a matrix with npop rows and ngen columns, with every entry = init
  for(i in 1:(ngen-1))
  {
    freq[,i+1] = rbinom(npop,nall,mu+(1-2*mu)*freq[,i]/nall) # create the 2,3,4, ... ,ngen generations, in each of the npop populations, by random sampling from the binomial distribution.  The expected allele fraction in the new generation is the allele fraction in the current generation (freq[,i]/nall) modified by the effect of mutation
  }
  par(mar=c(4,4.5,2.5,0.5));matplot(t(freq),type="l",lty=1,lwd=2,xlab="generations",ylab="allele count",main=paste("WF model in", npop," populations, mu = ",mu))
  return(freq[,ngen])
}

# het function
het.calc=function(x)
{
  prop.0=sum(x==0)/length(x)
  prop.1=sum(x==1)/length(x)
  H=1-prop.0^2-prop.1^2
  return(H)
}

hwe.calc=function(x)
{
  geno.X=apply(cbind(x[seq(1,length(x),2)],x[seq(2,length(x),2)]),1,sum)
  p_X=sum(x==0)/length(x)
  expected.vals=c(p_X^2,2*p_X*(1-p_X),(1-p_X)^2)*length(geno.X)
  observed.vals=c(sum(geno.X==0),sum(geno.X==1),sum(geno.X==2))
  pearson.chi.sq=sum(((observed.vals-expected.vals)^2)/expected.vals)
  p.val=1-pchisq(pearson.chi.sq,1)
  #return(pearson.chi.sq)
  return(p.val)
}


##Question 1

#(a) Imagine a population evolves for 500 generations under the Wright-Fisher model, with a starting allele frequency of 20%. What would the size of that population have to be such that roughly half of all SNPs with this starting frequency were either fixed or lost?

#By changing parametres, when 755<nall<760, there would be half of all SNPs fixed/lost.

#ac=wf(npop =10000, ngen=500,init=200, nall = 1000); table(ac)
#ac=wf(npop =10000, ngen=500,init=100, nall = 500); table(ac)
#ac=wf(npop =10000, ngen=500,init=140, nall = 700); table(ac)
#ac=wf(npop =10000, ngen=500,init=150, nall = 750); table(ac)
#ac=wf(npop =10000, ngen=500,init=154, nall = 770); table(ac)
#ac=wf(npop =10000, ngen=500,init=152, nall = 760); table(ac)

ac=wf(npop =10000, ngen=500,init=151, nall = 755); table(ac)

ac=wf(npop =10000, ngen=500,init=170, nall = 850); table(ac) # why not this one???


#(b)Simulate three populations with the same starting frequency of 20% for 500 generations, but with each population having different population sizes: 1000, 2000, and 20000. Calculate Heterozygosity (H) for each population at the end of the 500 generations. What do you see?
ac1=wf(npop =10000, ngen=500,init=200, nall = 1000); 

#calculate heterozygosity
H.ac1<-het.calc(ac1)

ac2=wf(npop =10000, ngen=500,init=400, nall = 2000); #table(ac)

#calculate heterozygosity

H.ac2 <-het.calc(ac2)

ac3=wf(npop =10000, ngen=500,init=4000, nall = 20000); #table(ac)

#calculate heterozygosity

H.ac3 <-het.calc(ac3)

#calculated heterozygosity for the three population
H.ac1
H.ac2
H.ac3



#(c)Calculate FST between the 3 populations in (b). Repeat this, but simulating for different numbers of generations (though with the same number of generations for each population). What do you see, and what implications does this have for real data?

#50 generations
ac1=wf(npop =10000, ngen=50,init=200, nall = 1000); #table(ac1)
ac2=wf(npop =10000, ngen=50,init=400, nall = 2000); #table(ac2)
ac3=wf(npop =10000, ngen=50,init=4000, nall = 20000); #table(ac3)


p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst50<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst50

#100 generations
ac1=wf(npop =10000, ngen=100,init=200, nall = 1000); #table(ac1)
ac2=wf(npop =10000, ngen=100,init=400, nall = 2000); #table(ac2)
ac3=wf(npop =10000, ngen=100,init=4000, nall = 20000); #table(ac3)


p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst100<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst100



#250 generations
ac1=wf(npop =10000, ngen=250,init=200, nall = 1000); #table(ac1)
ac2=wf(npop =10000, ngen=250,init=400, nall = 2000); #table(ac2)
ac3=wf(npop =10000, ngen=250,init=4000, nall = 20000); #table(ac3)


p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst250<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst250


#500 generations
ac1=wf(npop =10000, ngen=500,init=200, nall = 1000); #table(ac1)
ac2=wf(npop =10000, ngen=500,init=400, nall = 2000); #table(ac2)
ac3=wf(npop =10000, ngen=500,init=4000, nall = 20000); #table(ac3)
p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst500<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst500



#1000 generations
ac1=wf(npop =10000, ngen=1000,init=200, nall = 1000); #table(ac1)
ac2=wf(npop =10000, ngen=1000,init=400, nall = 2000); #table(ac2)
ac3=wf(npop =10000, ngen=1000,init=4000, nall = 20000); #table(ac3)
p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst1000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst1000


#2000 generations
ac1=wf(npop =10000, ngen=2000,init=200, nall = 1000); 
ac2=wf(npop =10000, ngen=2000,init=400, nall = 2000); #table(ac)
ac3=wf(npop =10000, ngen=2000,init=4000, nall = 20000); #table(ac)
p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst2000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst2000

#3000 generations
ac1=wf(npop =10000, ngen=3000,init=200, nall = 1000); 
ac2=wf(npop =10000, ngen=3000,init=400, nall = 2000); #table(ac)
ac3=wf(npop =10000, ngen=3000,init=4000, nall = 20000); #table(ac)
p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst3000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst3000

#6000 generations
ac1=wf(npop =10000, ngen=6000,init=200, nall = 1000); 
ac2=wf(npop =10000, ngen=6000,init=400, nall = 2000); #table(ac)
ac3=wf(npop =10000, ngen=6000,init=4000, nall = 20000); #table(ac)
p<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3
p

Fst6000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst6000



Fst250
Fst500
Fst1000
Fst2000
Fst3000



#(d)Now repeat (b) using simulations with selection. What do you see?

selectSIM <- function(ngen=500,nall=100,init=1,s=0)
{
  ## note: "init" = number of copies of allele at generation 0
  ## note: "nall" = number of diploids (not haploids)
  p = init/(2*nall)
  pop.samp = sample(c(2,1,0),nall,prob=c(p^2,2*p*(1-p),(1-p)^2),replace=TRUE)
  freq = matrix(NA,4,ngen) ## first 3 rows = genotype counts; 4th row = allele counts
  freq[1:3,1] = c(sum(pop.samp==2),sum(pop.samp==1),sum(pop.samp==0))
  freq[4,1] = 2*freq[1,1]+freq[2,1]
  for(i in 1:(ngen-1))
  {
    p = freq[4,i]/(2*(freq[1,i]+freq[2,i]+freq[3,i]))
    prob.vals = c((1+s)*p^2,2*p*(1-p),(1-p)^2)
    prob.vals[prob.vals<0]=0
    c.val = sum(prob.vals)
    prob.vals = prob.vals/c.val
    pop.samp = sample(c(2,1,0),nall,prob=prob.vals,replace=TRUE)
    freq[1:3,i+1] = c(sum(pop.samp==2),sum(pop.samp==1),sum(pop.samp==0))
    freq[4,i+1] = 2*freq[1,i+1]+freq[2,i+1]
  }
  to.plot='no'
  if (to.plot=='yes')
  {
    ylim.plot=c(0,2*nall)
    par(mar=c(4,4.5,2.5,0.5));matplot(t(freq),type="l",col=1:4,lty=1,lwd=5,xlab="generations",ylab="genotype/allele count",main=paste("init = ",init,", nall = ",nall,", s = ",s),ylim=ylim.plot,cex.main=2,cex.lab=1.5,cex.axis=1)
    abline(h=nall,lty=3)
    legend(0,max(ylim.plot),legend=c("AA","Aa","aa","A"),col=1:4,lwd=4,cex=1.2)
  }
  return(freq[,ngen])
}




#calculate heterozygosity
select1 <-selectSIM(ngen=500,nall=1000,init=200,s = 0.01)


#calculate heterozygosity
select2 <-selectSIM(ngen=500,nall=2000,init=400,s = 0.01)


#calculate heterozygosity
select3 <-selectSIM(ngen=500,nall=20000,init=4000,s = 0.01)


select1[1]


het_select1=1-((select1[4]/2000)^2)-(1-(select1[4]/2000))^2
het_select2=1-((select2[4]/4000)^2)-(1-(select2[4]/4000))^2
het_select3=1-((select3[4]/40000)^2)-(1-(select3[4]/40000))^2

het_select1
het_select2
het_select3



p<-select1[4]/2000
q<-1-(select1[4]/2000)
twopq<-2*p*q
#p
#q
twopq
h1<-1-(p^2)-(q^2)
h1



p<-select2[4]/4000
q<-1-(select2[4]/4000)
twopq<-2*p*q
#p
#q
twopq
h2<-1-(p^2)-(q^2)
h2



p<-select3[4]/40000
q<-1-(select3[4]/40000)
twopq<-2*p*q
#p
#q
twopq
h3<-1-(p^2)-(q^2)
h3



#(e)Now repeat (c), but adding selection to only one population. What do you see, and what implications does this have for real data?

ac1=selectSIM(ngen=50,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=50,nall=2000,init=400,s = 0.01); #table(ac)
ac3=selectSIM(ngen=50,nall=20000,init=4000,s = 0.01); #table(ac)
p50<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst50<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst50

#100 generations
ac1=selectSIM(ngen=100,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=100,nall=2000,init=400); #table(ac)
ac3=selectSIM(ngen=100,nall=20000,init=4000); #table(ac)
p100<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst100<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst100

#250 generations
ac1=selectSIM(ngen=250,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=250,nall=2000,init=400); #table(ac)
ac3=selectSIM(ngen=250,nall=20000,init=4000); #table(ac)
p250<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst250<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst250


#500 generations
ac1=selectSIM(ngen=500,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=500,nall=2000,init=400); #table(ac)
ac3<-selectSIM(ngen=500,nall=20000,init=4000); #table(ac)
p500<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst500<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst500



#1000 generations
ac1=selectSIM(ngen=1000,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=1000,nall=2000,init=400); #table(ac)
ac3<-selectSIM(ngen=1000,nall=20000,init=4000); #table(ac)
p1000<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst1000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst1000

#2000 generations
ac1=selectSIM(ngen=2000,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=2000,nall=2000,init=400); #table(ac)
ac3<-selectSIM(ngen=2000,nall=20000,init=4000); #table(ac)
p2000<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst2000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst2000


#3000 generations
ac1=selectSIM(ngen=3000,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=3000,nall=2000,init=400); #table(ac)
ac3<-selectSIM(ngen=3000,nall=20000,init=4000); #table(ac)
p3000<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst3000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

Fst3000

#6000 generations
ac1=selectSIM(ngen=6000,nall=1000,init=200,s = 0.01); 
ac2=selectSIM(ngen=6000,nall=2000,init=400); #table(ac)
ac3<-selectSIM(ngen=6000,nall=20000,init=4000); #table(ac)
p6000<-(mean(ac1)/1000+mean(ac2)/2000+mean(ac3)/20000)/3


Fst6000<-(((mean(ac1)/1000-p)^2)+((mean(ac2)/2000-p)^2)+((mean(ac3)/20000-p)^2))/(3*p*(1-p))

p50
p100
p250
p500
p1000
p2000
p3000
p6000



Fst50
Fst100
Fst250
Fst500
Fst1000
Fst2000
Fst3000
Fst6000






#Question 2
het.calc=function(x)
{
  prop.0=sum(x==0)/length(x)
  prop.1=sum(x==1)/length(x)
  H=1-prop.0^2-prop.1^2
  return(H)
}

hwe.calc=function(x)
{
  geno.X=apply(cbind(x[seq(1,length(x),2)],x[seq(2,length(x),2)]),1,sum)
  p_X=sum(x==0)/length(x)
  expected.vals=c(p_X^2,2*p_X*(1-p_X),(1-p_X)^2)*length(geno.X)
  observed.vals=c(sum(geno.X==0),sum(geno.X==1),sum(geno.X==2))
  pearson.chi.sq=sum(((observed.vals-expected.vals)^2)/expected.vals)
  p.val=1-pchisq(pearson.chi.sq,1)
  #return(pearson.chi.sq)
  return(p.val)
}

#---------------------------------------------------------------------
# A
ceu = t(read.table(file.choose()))
chb = t(read.table(file.choose()))
gih = t(read.table(file.choose()))
jpt = t(read.table(file.choose()))
lwk = t(read.table(file.choose()))
yri = t(read.table(file.choose()))

table(ceu)
table(chb)
table(gih)
table(jpt)
table(lwk)
table(yri)

allele_frequency.ceu = table(ceu)[2]/(table(ceu)[1]+table(ceu)[2])
allele_frequency.chb = table(chb)[2]/(table(chb)[1]+table(chb)[2])
allele_frequency.gih = table(gih)[2]/(table(gih)[1]+table(gih)[2])
allele_frequency.jpt = table(jpt)[2]/(table(jpt)[1]+table(jpt)[2])
allele_frequency.lwk = table(lwk)[2]/(table(lwk)[1]+table(lwk)[2])
allele_frequency.yri = table(yri)[2]/(table(yri)[1]+table(yri)[2])

allele_frequency.ceu
allele_frequency.chb
allele_frequency.gih
allele_frequency.jpt
allele_frequency.lwk
allele_frequency.yri

# E (a)
ceu2 = t(read.table(file.choose()))
chb2 = t(read.table(file.choose()))
gih2 = t(read.table(file.choose()))
jpt2 = t(read.table(file.choose()))
lwk2 = t(read.table(file.choose()))
yri2 = t(read.table(file.choose()))

table(ceu2)
table(chb2)
table(gih2)
table(jpt2)
table(lwk2)
table(yri2)

allele_frequency.ceu2 = table(ceu2)[2]/(table(ceu2)[1]+table(ceu2)[2])
allele_frequency.chb2 = table(chb2)[2]/(table(chb2)[1]+table(chb2)[2])
allele_frequency.gih2 = table(gih2)[2]/(table(gih2)[1]+table(gih2)[2])
allele_frequency.jpt2 = table(jpt2)[2]/(table(jpt2)[1]+table(jpt2)[2])
allele_frequency.lwk2 = table(lwk2)[2]/(table(lwk2)[1]+table(lwk2)[2])
allele_frequency.yri2 = table(yri2)[2]/(table(yri2)[1]+table(yri2)[2])

allele_frequency.ceu2
allele_frequency.chb2
allele_frequency.gih2
allele_frequency.jpt2
allele_frequency.lwk2
allele_frequency.yri2

# D
# Code ------------------------------------------
# Using the data x,y from any two SNPs
d.prime.calc=function(x,y)
{
  D.00=length(x[x==0 & y==0])/length(x)-(length(x[x==0])/
                                           length(x))*(length(y[y==0])/length(y))
  D.minus=min((length(x[x==1])/length(x))*(length(y[y==1])/length(y)),
              (length(x[x==0])/length(x))*(length(y[y==0])/length(y)))
  D.plus=min((length(x[x==1])/length(x))*(length(y[y==0])/length(y)),
             (length(x[x==0])/length(x))*(length(y[y==1])/length(y)))
  if (D.00>=0) D.prime=D.00/D.plus
  if (D.00<0) D.prime=D.00/D.minus
  return(abs(D.prime))
}
#-------------------------------------------------------------------------------
# ceu
num.snps=dim(ceu)[2] 
min.allelefreq.ceu=D.prime.ceu=cor.ceu=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.ceu[i,j]=d.prime.calc(ceu[,i],ceu[,j])
    cor.ceu[i,j]=(cor(ceu[,i],ceu[,j]))^2
    min.allelefreq.ceu[i,j]=min(c(sum(ceu[,i]==0),sum(ceu[,i]==1),
                              sum(ceu[,j]==0),sum(ceu[,j]==1))/dim(ceu)[1])
  } 
}

dim(min.allelefreq.ceu)
dim(D.prime.ceu)
dim(cor.ceu)
#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.ceu=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.ceu[i]=mean(D.prime.ceu[min.allelefreq.ceu>allelefreq.bins[i]
& min.allelefreq.ceu<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.ceu


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.ceu=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.ceu[i]=mean(cor.ceu[min.allelefreq.ceu>allelefreq.bins[i]
& min.allelefreq.ceu<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.ceu
#--------------------------------------------------------------------------------
# chb
num.snps=dim(chb)[2] 
min.allelefreq.chb=D.prime.chb=cor.chb=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.chb[i,j]=d.prime.calc(chb[,i],chb[,j])
    cor.chb[i,j]=(cor(chb[,i],chb[,j]))^2
    min.allelefreq.chb[i,j]=min(c(sum(chb[,i]==0),sum(chb[,i]==1),
                                  sum(chb[,j]==0),sum(chb[,j]==1))/dim(chb)[1])
  } 
}

dim(min.allelefreq.chb)
dim(D.prime.chb)
dim(cor.chb)
#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.chb=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.chb[i]=mean(D.prime.chb[min.allelefreq.chb>allelefreq.bins[i]
                                       & min.allelefreq.chb<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.chb


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.chb=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.chb[i]=mean(cor.chb[min.allelefreq.chb>allelefreq.bins[i]
                               & min.allelefreq.chb<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.chb
#--------------------------------------------------------------------------------
# gih
num.snps=dim(gih)[2] 
min.allelefreq.gih=D.prime.gih=cor.gih=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.gih[i,j]=d.prime.calc(gih[,i],gih[,j])
    cor.gih[i,j]=(cor(gih[,i],gih[,j]))^2
    min.allelefreq.gih[i,j]=min(c(sum(gih[,i]==0),sum(gih[,i]==1),
                                  sum(gih[,j]==0),sum(gih[,j]==1))/dim(gih)[1])
  } 
}

dim(min.allelefreq.gih)
dim(D.prime.gih)
dim(cor.gih)


#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.gih=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.gih[i]=mean(D.prime.gih[min.allelefreq.gih>allelefreq.bins[i]
                                       & min.allelefreq.gih<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.gih


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.gih=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.gih[i]=mean(cor.gih[min.allelefreq.gih>allelefreq.bins[i]
                               & min.allelefreq.gih<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.gih
#--------------------------------------------------------------------------------
# jpt
num.snps=dim(jpt)[2] 
min.allelefreq.jpt=D.prime.jpt=cor.jpt=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.jpt[i,j]=d.prime.calc(jpt[,i],jpt[,j])
    cor.jpt[i,j]=(cor(jpt[,i],jpt[,j]))^2
    min.allelefreq.jpt[i,j]=min(c(sum(jpt[,i]==0),sum(jpt[,i]==1),
                                  sum(jpt[,j]==0),sum(jpt[,j]==1))/dim(jpt)[1])
  } 
}

dim(min.allelefreq.jpt)
dim(D.prime.jpt)
dim(cor.jpt)
#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.jpt=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.jpt[i]=mean(D.prime.jpt[min.allelefreq.jpt>allelefreq.bins[i]
                                       & min.allelefreq.jpt<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.jpt


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.jpt=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.jpt[i]=mean(cor.jpt[min.allelefreq.jpt>allelefreq.bins[i]
                               & min.allelefreq.jpt<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.jpt
#--------------------------------------------------------------------------------
# lwk
num.snps=dim(lwk)[2] 
min.allelefreq.lwk=D.prime.lwk=cor.lwk=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.lwk[i,j]=d.prime.calc(lwk[,i],lwk[,j])
    cor.lwk[i,j]=(cor(lwk[,i],lwk[,j]))^2
    min.allelefreq.lwk[i,j]=min(c(sum(lwk[,i]==0),sum(lwk[,i]==1),
                                  sum(lwk[,j]==0),sum(lwk[,j]==1))/dim(lwk)[1])
  } 
}

dim(min.allelefreq.lwk)
dim(D.prime.lwk)
dim(cor.lwk)
#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.lwk=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.lwk[i]=mean(D.prime.lwk[min.allelefreq.lwk>allelefreq.bins[i]
                                       & min.allelefreq.lwk<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.lwk


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.lwk=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.lwk[i]=mean(cor.lwk[min.allelefreq.lwk>allelefreq.bins[i]
                               & min.allelefreq.lwk<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.lwk
#--------------------------------------------------------------------------------
# yri
num.snps=dim(yri)[2] 
min.allelefreq.yri=D.prime.yri=cor.yri=matrix(NA,nrow=num.snps,ncol=num.snps) 
for (i in 1:(num.snps-1))
{
  for(j in (i+1):num.snps)
  {
    D.prime.yri[i,j]=d.prime.calc(yri[,i],yri[,j])
    cor.yri[i,j]=(cor(yri[,i],yri[,j]))^2
    min.allelefreq.yri[i,j]=min(c(sum(yri[,i]==0),sum(yri[,i]==1),
                                  sum(yri[,j]==0),sum(yri[,j]==1))/dim(yri)[1])
  } 
}

dim(min.allelefreq.yri)
dim(D.prime.yri)
dim(cor.yri)
#Average
allelefreq.bins=seq(0,0.5,by=0.01) 
mean.D.prime.yri=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.D.prime.yri[i]=mean(D.prime.yri[min.allelefreq.yri>allelefreq.bins[i]
                                       & min.allelefreq.yri<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.D.prime.yri


allelefreq.bins=seq(0,0.5,by=0.01) 
mean.cor.yri=rep(NA,length(allelefreq.bins)-1) 
for (i in 1:(length(allelefreq.bins)-1))
{
  mean.cor.yri[i]=mean(cor.yri[min.allelefreq.yri>allelefreq.bins[i]
                               & min.allelefreq.yri<allelefreq.bins[i+1]],na.rm=TRUE) 
}
mean.cor.yri

#Plot-------------------------------------------------------------------------------------
df.plot.ceu <- data.frame("D.prime"= mean.D.prime.ceu,"R_square"= mean.cor.ceu, "MAF"=seq(0.01,0.5,by=0.01))
df.plot.chb <- data.frame("D.prime"= mean.D.prime.chb,"R_square"= mean.cor.chb, "MAF"=seq(0.01,0.5,by=0.01))
df.plot.gih <- data.frame("D.prime"= mean.D.prime.gih,"R_square"= mean.cor.gih, "MAF"=seq(0.01,0.5,by=0.01))
df.plot.jpt <- data.frame("D.prime"= mean.D.prime.jpt,"R_square"= mean.cor.jpt, "MAF"=seq(0.01,0.5,by=0.01))
df.plot.lwk <- data.frame("D.prime"= mean.D.prime.lwk,"R_square"= mean.cor.lwk, "MAF"=seq(0.01,0.5,by=0.01))
df.plot.yri <- data.frame("D.prime"= mean.D.prime.yri,"R_square"= mean.cor.yri, "MAF"=seq(0.01,0.5,by=0.01))

library(ggplot2)
ggplot(df.plot.ceu, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("CEU") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.chb, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("CHB") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.gih, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("GIH") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.jpt, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("JPT") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.lwk, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("LWK") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.yri, aes(x=MAF, y=D.prime)) + geom_point(shape = 17, color="steelblue3") + ylim(0, 1) + ggtitle("YRI") +  theme(plot.title = element_text(hjust = 0.5))

ggplot(df.plot.ceu, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("CEU") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.chb, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("CHB") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.gih, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("GIH") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.jpt, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("JPT") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.lwk, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("LWK") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(df.plot.yri, aes(x=MAF, y=R_square)) + geom_point(shape = 17, color="pink2") + ylim(0, 0.12) + ggtitle("YRI") +  theme(plot.title = element_text(hjust = 0.5))






# Excercise 2
CEU = read.table(file.choose())
JPT_CHB = read.table(file.choose())
YRI = read.table(file.choose())
#----------------------------------
#Heterozygosity-------------------------------
het.calc(CEU[,1])

num.snps=dim(CEU)[2]
het.CEU=rep(NA,num.snps)
for (i in 1:num.snps) het.CEU[i]=het.calc(CEU[,i])
het.CEU
summary(het.CEU)

num.snps=dim(JPT_CHB)[2]
het.JPT_CHB=rep(NA,num.snps)
for (i in 1:num.snps) het.JPT_CHB[i]=het.calc(JPT_CHB[,i])
het.JPT_CHB
summary(het.JPT_CHB)

num.snps=dim(YRI)[2]
het.YRI=rep(NA,num.snps)
for (i in 1:num.snps) het.YRI[i]=het.calc(YRI[,i])
het.YRI
summary(het.YRI)

#Hardy-Weinburg Equilibrium (HWE)--------------------------
hwe.calc(CEU[,1])
all=rbind(CEU,JPT_CHB,YRI)

hwe.snps.CEU=hwe.snps.JPT_CHB=hwe.snps.YRI=hwe.snps.all=rep(NA,num.snps)
all=rbind(CEU,JPT_CHB,YRI)
for (i in 1:num.snps) hwe.snps.CEU[i]=hwe.calc(CEU[,i])
for (i in 1:num.snps) hwe.snps.JPT_CHB[i]=hwe.calc(JPT_CHB[,i])
for (i in 1:num.snps) hwe.snps.YRI[i]=hwe.calc(YRI[,i])
for (i in 1:num.snps) hwe.snps.all[i]=hwe.calc(all[,i])

par(mfrow=c(2,2)) ## plots 2 figures per row, and 2 per column
hist(hwe.snps.CEU,main="CEU",breaks=seq(0,1,0.05),xlab="HWE p-values
   across SNPs")
hist(hwe.snps.JPT_CHB,main="JPT_CHB",breaks=seq(0,1,0.05),xlab="HWE p-values
   across SNPs")
hist(hwe.snps.YRI,main="YRI",breaks=seq(0,1,0.05),xlab="HWE p-values
   across SNPs")
hist(hwe.snps.all,main="ALL",breaks=seq(0,1,0.05),xlab="HWE p-values
   across SNPs")

par(mfrow=c(1,1)) ## goes back to one figure total
boxplot(hwe.snps.CEU,hwe.snps.JPT_CHB,hwe.snps.YRI,hwe.snps.all,
        ylab="HWE p-values",names=c("CEU","JPT","YRI","ALL"))










