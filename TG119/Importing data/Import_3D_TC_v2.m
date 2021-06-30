
%% Rotina para  importação de imagens DICOM, e compilação em estrutura 3D

clear all
close all

[filename, pathname] = uigetfile('*','Selecione as imagens Tomograficas','MultiSelect','on');
num= length(filename);
cont=1;

for aaaa = 1:num
    close
    xinfo(aaaa)=dicominfo([pathname,char(filename(aaaa))]);
    x=dicomread([pathname,char(filename(aaaa))]);
    
    TC_full(:,:,cont)=x; %Volume com os slices agrupados
    cont=cont+1;
end

%%
%Etapa para leitura das estruturas 
%dados estão no cabeçalho da imagem dicom, por isso a função dicominfo.

[filename1, pathname1] = uigetfile('*','Selecione o conjunto de estruturas');
estrutura=dicominfo([pathname1,filename1]);

contornos=readRTstructures(estrutura,xinfo);

%% Estruturas Prostata

Body=contornos(1).Segmentation;
PTV=contornos(2).Segmentation;
Prostate=contornos(3).Segmentation;
Urinary_bladder=contornos(4).Segmentation;
Rectum=contornos(5).Segmentation;


figure, isosurface(Body)
hold on
alpha 0.3
isosurface(PTV)
alpha .7
isosurface(Prostate)

isosurface(Rectum)
isosurface(Urinary_bladder)