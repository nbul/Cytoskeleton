clear variables
clc
close all

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

Data = readtable('Summary_pulled.csv');
Data.directionActin
Data.Direction_cell
Data.MTDirection

NormMT = Data.MTDirection - Data.Direction_cell;
NormMT(NormMT<-90) = 180 + NormMT(NormMT<-90);
NormMT(NormMT>90) = -180 + NormMT(NormMT>90);
NormA = Data.directionActin - Data.Direction_cell;
NormA(NormA<-90) = 180 + NormA(NormA<-90);
NormA(NormA>90) = -180 + NormA(NormA>90);

NormA2 = NormA;
NormMT2 = NormMT;
NormA2(NormA<0) = -NormA2(NormA<0);
NormMT2(NormA<0) = -NormMT2(NormA<0);
 %180+NormMT2((NormMT2-NormA2)<-90);

NormMT215 = NormMT2;
NormMT215(abs((NormMT2-NormA2))>90 | NormA2<15) =[];

NormMT240 = NormMT2;
NormMT240(abs((NormMT2-NormA2))>90 | NormA2<25) = [];

NormMT2(abs((NormMT2-NormA2))>90) =[];
cd(currdir);