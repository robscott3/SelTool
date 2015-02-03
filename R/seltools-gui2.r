
library(tcltk)
library(FLCore)

## GLOBAL VARIABLES
memcount <- 0
res      <- popsel()$gearsel

## FUNCTIONS

popsel <- function(fmult=c(0.4,0.2,0.1), rcts=rep(1000,3), mix=0.05, symmetric=T){

  # define the dimensions
  ages  <- 1:15
  areas <- 1:length(fmult)

  # define static parameters
  nm    <- 0.2
  a50   <- 4
  b     <- 1.5
  sel   <- FLQuant(1/(1+exp(-b*(ages-a50))),dimnames=list(age=ages))

  # movement matrix
  mixmat <- array(mix, dim=c(length(areas),length(areas)))
  if(symmetric)
    diag(mixmat) <- 1-(mix*(length(areas)-1))

  if(!symmetric){
    mixmat <- diag(c(rep(1-mix,length(fmult)-1),1))
    mixmat[lower.tri(mixmat)][cumsum(1:(length(fmult)-1))] <- mix
  }
  # calculate N at age
  Ns <- FLQuant(NA, dimnames=list(age=ac(ages), area=ac(areas)))
  Ns[1,] <- rcts
  for(j in ages[-min(ages)]) {
    for(k in areas) {
      Ns[j,,,,k,] <-sum(Ns[j-1,,,,,] * mixmat[k,] * exp(-nm -fmult*c(sel[j-1])))
    }
  }

  # population selection
  popF   <- log(areaSums(Ns)[-max(ages)]/areaSums(Ns)[-min(ages)])-nm
  popsel <- popF/max(popF)
  return(FLQuants(popsel=popsel, gearsel=sel[-max(ages)]))
}
  

OnClear <- function() {
  mem  <- as.numeric(tclvalue(MEM))
  memcount <<- 2
  res <<- propagate(FLQuant(NA, dimnames=dimnames(popsel()$gearsel)),mem+1)
  iter(res,1) <<- popsel()$gearsel
  plot(xyplot(data~age, groups=iter, data=res, type='l', col='black', lty=2))
}

OnAdd  <- function() {
  fs <- c(as.numeric(tclvalue(F1)), as.numeric(tclvalue(F2)), as.numeric(tclvalue(F3)))
  rs <- c(as.numeric(tclvalue(R1)), as.numeric(tclvalue(R2)), as.numeric(tclvalue(R3)))
  mx <- as.numeric(tclvalue(MIX))
  sym <- as.numeric(tclvalue(rbvalue))==1
  
  iter(res, memcount) <<- popsel(fs, rs, mx, sym)$popsel
  plot(xyplot(data~age, groups=iter, data=res, type='l', 
       col=c('black',rep('grey',memcount-2),'black'), lty=c(2,rep(1,memcount-1))))

  memcount <<- memcount + 1
}

OnReset <- function(){
  tclvalue(F1) <- '0.4'
}

## POINTY CLICKY STUFF

tt <- tktoplevel()
tkwm.title(tt, text='population selection tool')

F1 <- tclVar(0.4); F2 <- tclVar(0.2); F3 <- tclVar(0.1)
R1 <- tclVar(1000); R2 <- tclVar(1000); R3 <- tclVar(1000)
MIX <- tclVar(0.05)
MEM <- tclVar(25)
rbvalue <- tclVar(TRUE)


f1 <- tkentry(tt, width=6, textvariable=F1)
f2 <- tkentry(tt, width=6, textvariable=F2)
f3 <- tkentry(tt, width=6, textvariable=F3)
r1 <- tkentry(tt, width=6, textvariable=R1)
r2 <- tkentry(tt, width=6, textvariable=R2)
r3 <- tkentry(tt, width=6, textvariable=R3)
mix <- tkentry(tt, width=6, textvariable=MIX)
mem <- tkentry(tt, width=6, textvariable=MEM)
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)

tkconfigure(rb1, variable=rbvalue, value=TRUE)
tkconfigure(rb2, variable=rbvalue, value=FALSE)

addbutt    <- tkbutton(tt, text='Add',    command=OnAdd)
clearbutt  <- tkbutton(tt, text='Clear',  command=OnClear)
cancelbutt <- tkbutton(tt, text='Cancel', command=function()tkdestroy(tt))

tkgrid(tklabel(tt, width=15, text=''))
tkgrid(tklabel(tt, width=15, text='Areas'), tklabel(tt, width=5, text='Area 1'),
                                            tklabel(tt, width=5, text='Area 2'),
                                            tklabel(tt, width=5, text='Area 3'),
                                            tklabel(tt, width=5, text=''))
tkgrid(tklabel(tt, width=15, text='F'), f1, f2, f3)
tkgrid(tklabel(tt, width=15, text='Recruitment'), r1, r2, r3) 
tkgrid(tklabel(tt, width=15, text='Mixing'), mix)
tkgrid(tklabel(tt, width=15, text='Symmetric'), tklabel(tt, width=5, text='True'), rb1)
tkgrid(tklabel(tt, width=15, text=''),          tklabel(tt, width=5, text='False'), rb2)
tkgrid(tklabel(tt, width=15, text=''))
tkgrid(tklabel(tt, width=15, text='Memory'), mem)
tkgrid(tklabel(tt, width=15, text=''), addbutt, clearbutt, cancelbutt)
tkgrid(tklabel(tt, width=15, text=''))

