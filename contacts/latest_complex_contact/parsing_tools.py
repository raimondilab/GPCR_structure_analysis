#!/usr/bin/env python3


import os, sys, operator, math, gzip

################################################################################################################################
####### A collection of python functions to parse various formats of text files ######
#################################Written by Francesco Raimondi##################################################################
################################################################################################################################

f="\t"


def GetID(aa):
  labels=""
  resid=""
  flag=0
  flag2=0
  cc=0
  for c in aa:
    #print resid
    cc+=1
    if c.isdigit() == False:
      labels+=c
      flag2+=1
    elif c.isdigit() == True:
      flag=1
      resid+=c
    if c.isdigit() == False and flag==1 :
      #print resid
      if int(resid)>0:
        return resid
        break
    elif cc==len(aa) and flag==1 :
      #print resid
      if int(resid)>0:
        return resid
        break
    

