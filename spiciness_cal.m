%% INPUT
clear
addpath ./gsw/
addpath ./gsw/library
addpath ./gsw/thermodynamics_from_t
load('./input_index.mat');
load('./salt_temp.mat');
%%
depth(depth==10)=nan;
for ii=1:size(Cs_r,1)
    waterdepth(:,:,ii)=depth.*Cs_r(ii);
end
for ii=1:size(waterdepth,3)
    for i=1:size(waterdepth,1)
        for j=1:size(waterdepth,2)
            if (waterdepth(i,j,ii)<0)&&(mask_rho(i,j)~=0)
                valid_num(i,j,ii)=1;
            else
                valid_num(i,j,ii)=0;
            end
        end
    end
end
valid=valid_num~=0;
sr = salt;
sr(valid)= gsw_SR_from_SP(salt(valid)); 
sr(~valid) = NaN;
[npi,npj,npk]=size(salt);
ct = zeros(npi,npj,npk);
ct(valid) = gsw_CT_from_pt(sr(valid),temp(valid));
ct(~valid) = NaN;
sigma_ref = gsw_rho_CT_exact(sr,ct,0)-1000.; 
sigma_ref(~valid) = NaN; 
%% Compute Huang's spicity variable
pi_ref = zeros(npi,npj,npk);
pi_ref(:) = NaN; 
phuang = [0.,500.,1000.,2000.,3000.,4000.,5000.];
for k=1:npk
    disp(k) 
    for j=1:npj
        for i=1:npi
            if sr(i,j,k)>10.&&sr(i,j,k)<40.&&ct(i,j,k)>-2.&&ct(i,j,k)<40.
                pi_ref(i,j,k)=gsw_pspi(sr(i,j,k),ct(i,j,k),phuang(1));
            end
        end
    end
end
%%
load('input_index2.mat');
h_max=max(max(h));
[aa,bb]=size(h);
zz=[-9000:500:-500,-450:50:-5];
saltz(1:aa,1:bb,1:size(zz,2))=nan;
tempz(1:aa,1:bb,1:size(zz,2))=nan;
piz(1:aa,1:bb,1:size(zz,2))=nan;
sigmaz(1:aa,1:bb,1:size(zz,2))=nan;
for i=1:aa
    for j=1:bb
           z_r=((h(i,j)+zeta(i,j)).*Cs_r(1:30))+zeta(i,j);
           saltz(i,j,:)=interp1(z_r,squeeze(salt(i,j,:)),zz);
           tempz(i,j,:)=interp1(z_r,squeeze(temp(i,j,:)),zz);
           piz(i,j,:)=interp1(z_r,squeeze(pi_ref(i,j,:)),zz);%spiciness
           sigmaz(i,j,:)=interp1(z_r,squeeze(sigma_ref(i,j,:)),zz);%density
    end
end
save
