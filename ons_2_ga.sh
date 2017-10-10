#mkdir ons
#cd ons
##!/bin/bash

#outfile="diss_ons2.txt"
#rm $outfile
period=20
t_stop=1000
n_step=10000
d_ec=-0.0 #-100
d_ef=0.2 #-0.025 #-100
thr_rk=1e-7
tol_rk=1e-7
time_integr='rks'		#rks for rksuite, rkf for rkf45

scf_type='ga'
threshold_k=1e-8
threshold_kr=1e-8
threshold=1e-4
lpathnscf=1
loopr=1
nkx_plot=11 #801		#for G, G0 and arpes kr

T=0.01
alpha_min_k=0.1
lambda_par_k=0.1
alpha_mu_k=0.1
lk3d=0
alpha_max_k=${alpha_min_k}
lcheckcubic=0
nk_dos=${nkx_plot}
T_dos=${T_g}
iplane_xyz=3
nkx_chern=101
#outfile="size3.txt"
#sm2 =2x1 rec
#sm2b=1x1
#sm2c=2x1 non rec
label="ons_2_ga"

ec=0 #-6.74 #-1.46 #-2.4 #-2.465
ef=-0.1 #-0.576
tc=1 #1.9
V=0.6 #-0.3 #-0.5
tf=-0.1 #-0.02
nx=9
ny=${nx}
nx_kr=1
ny_kr=1
nkx_r=7 #`echo "63/${nx}"|bc -l`
nky_r=${nkx_r} #`echo "63/${ny}"|bc -l`
nz=1
lrec=0
d_e_rec=-100
d_ec_kr=0.0 #-0.2 #-0.05 #0.1 #0.1 #0.5 #0.05 #-100
d_ef_kr=${d_ec_kr} #-0.05 #0.1 #0.1 #0.5 #0.08 #-100
pbcz=0
ec0=1
nz_plot=3
nz_scf=1
phi_type="ons"
eta_7=1 #3 #0
eta2_f7=0 #-0.5 #-1 #-0.2
eta3_f7=0 #-0.5 #-0.5 #-0.1
eta1_d=1 #0.5
eta2_d=0 #-0.4 #-0.4 #-0.2
eta3_d=0 #-0.1
eta_v=1
eta_v2=0 #-0.1 #-0.2 #0.5
eta_v3=0

t_d1=1
t_f1=-0.19
lfullbz=0
stm_layer_ratio=0.1
alfa_nz=0.3

loopk=1
lscfk=1
loopkr=1
lscfkr=1
lscfloop_kr=1		#scf kr
lscfr=1
lscfloop=1
lnscfr=0
lgreen=0
lborn=0
lg0=0
lv_nscf=0		#uses nscf solution for V
lpathnscfr=0
lreadparam=0
lkshifted=0
nkx_r_nscf=10		#for arpes r
nky_r_nscf=10
lkshifted_r_nscf=0
lkred_r=1
ltr=1

lrhoktot=0
lbinning=1
lfixmu_k=0
lfixmu_kr=0
lfixmu=1
lmucf=1

T_g=0.003	#broadening for G
omega_min=-1
omega_max=1
delta_omega=0.01

nkpath_plot=100
lvacancy=0	#0
lvacancy4=0	#1
lvacancy3=0	#nz/2
lvacancy5=0	#nz-2
lvacancy2=0	#nz-1
ldisorder_r=0
ldisorder_r2=0
conc_vac=0.1
#T_max=1
T_max=`echo "${T}+0.00001"|bc -l`
#T_max=20
delta_t=0.01

lwfc=.false.
e_thr1=0.328
e_thr2=0.330

#other="_kr"

de=0

limpurity=1	#0
limpurity2=0	#1



b_init=0.48 #0.72
#b_init=0.7
#b_init=0.7
#lambda_init=`echo "${ef}-0.1"|bc -l`
lambda_init=-0.7 #-0.2
#lambda_init=`echo "${ef}-3.1"|bc -l`
#lambda_init=`echo "${ef}-6.1"|bc -l`
muinit=0.6 #-0.4 #-0.13


#############################################

#for ef in   -7 -6 -5 -4  -3 -2 -1.5 -1.4 -1.3 -1.2 -1.1  -1 -0.5 -0.2 -0.1
#for ef in  -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8
#for ef in    -1.4 -1.3 -1.2 -1.1  -1 -0.5 -0.2 -0.1 0.0 0.1 0.3 0.5 1 2 3 4
#for ef in      7 8
#for V in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10
#for V in 2.0 #2.1 #2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10.0 10.4 10.8 11.2 11.6 12.0 12.4 12.8 13.2 13.6 14.0
#for V in 11.2 11.6 12.0 12.4 12.8
#for V in 14.0 13.6 13.2 12.8 12.4 12.0 11.6 11.2 10.8 10.4 10.0 9.6 9.2 8.8 8.4 8.0 7.6 7.2 6.8 6.4 6.0 5.6 5.2 4.8 4.4 4.0 3.6 3.2 2.8 2.4 2.0 1.6 1.2 0.8 0.4
#for V in `seq 0.09 -0.01 0`
#for V in 1.6 #12 16
#for V in 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10
#for V in    #3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10
#for V in 6 #12 16
#for tf in   0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 # -0.1  # -7 -6 -5 -4  -3 -2 -1.5 -1.4 -1.3 -1.2 -1.1  -1 -0.5 -0.2 -0.1
#for de in 0 1 2 3 4 5
#do
#do
#do
#for d_ef in 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 -0.05 -0.1 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5 -0.55 -0.6 -0.65 -0.7 -0.75 -0.8 -0.85
for d_ef in 0.35
do
##############################################

label="ons2c_de${d_ef}"

ef0=`echo "${ef}+(${de})"|bc -l`


test ${lrhoktot} = 0 && lrhoktotb=.false.
test ${lrhoktot} = 1 && lrhoktotb=.true.
test ${ldisorder_r} = 0 && ldisorder_rb=.false.
test ${ldisorder_r} = 1 && ldisorder_rb=.true.
test ${ldisorder_r2} = 0 && ldisorder_r2b=.false.
test ${ldisorder_r2} = 1 && ldisorder_r2b=.true.
test ${loopk} = 0 && loopkb=.false.
test ${loopk} = 1 && loopkb=.true.
test ${lscfk} = 0 && lscfkb=.false.
test ${lscfk} = 1 && lscfkb=.true.
test ${lscfkr} = 0 && lscfkrb=.false.
test ${lscfkr} = 1 && lscfkrb=.true.
test ${loopkr} = 0 && loopkrb=.false.
test ${loopkr} = 1 && loopkrb=.true.
test ${lgreen} = 0 && lgreenb=.false.
test ${lgreen} = 1 && lgreenb=.true.
test ${lv_nscf} = 0 && lv_nscfb=.false.
test ${lv_nscf} = 1 && lv_nscfb=.true.
test ${lborn} = 0 && lbornb=.false.
test ${lborn} = 1 && lbornb=.true.
test ${loopr} = 0 && looprb=.false.
test ${loopr} = 1 && looprb=.true.
test ${lscfr} = 0 && lscfrb=.false.
test ${lscfr} = 1 && lscfrb=.true.
test ${lscfloop} = 0 && lscfloopb=.false.
test ${lscfloop} = 1 && lscfloopb=.true.
test ${lscfloop_kr} = 0 && lscfloop_krb=.false.
test ${lscfloop_kr} = 1 && lscfloop_krb=.true.
test ${lnscfr} = 0 && lnscfrb=.false.
test ${lnscfr} = 1 && lnscfrb=.true.
test ${lpathnscfr} = 0 && lpathnscfrb=.false.
test ${lpathnscfr} = 1 && lpathnscfrb=.true.
test ${lpathnscf} = 0 && lpathnscfb=.false.
test ${lpathnscf} = 1 && lpathnscfb=.true.
test ${lbinning} = 0 && lbinningb=.false.
test ${lbinning} = 1 && lbinningb=.true.
test ${lfixmu_kr} = 0 && lfixmu_krb=.false.
test ${lfixmu_kr} = 1 && lfixmu_krb=.true.
test ${lfixmu_k} = 0 && lfixmu_kb=.false.
test ${lfixmu_k} = 1 && lfixmu_kb=.true.
test ${lfixmu} = 0 && lfixmu_b=.false.
test ${lfixmu} = 1 && lfixmu_b=.true.
test ${lmucf} = 0 && lmucfb=.false.
test ${lmucf} = 1 && lmucfb=.true.
test ${lreadparam} = 0 && lreadparamb=.false.
test ${lreadparam} = 1 && lreadparamb=.true.
test ${lkshifted_r_nscf} = 0 && lkshifted_r_nscfb=.false.
test ${lkshifted_r_nscf} = 1 && lkshifted_r_nscfb=.true.
test ${lkshifted} = 0 && lkshiftedb=.false.
test ${lkshifted} = 1 && lkshiftedb=.true.
test ${lkred_r} = 0 && lkred_rb=.false.
test ${lkred_r} = 1 && lkred_rb=.true.
test ${ltr} = 0 && ltrb=.false.
test ${ltr} = 1 && ltrb=.true.
test ${pbcz} = 0 && pbczb=.false.
test ${pbcz} = 1 && pbczb=.true.
test ${lg0} = 0 && lg0b=.false.
test ${lg0} = 1 && lg0b=.true.
test ${lfullbz} = 0 && lfullbzb=.false.
test ${lfullbz} = 1 && lfullbzb=.true.
test ${lvacancy} = 0 && lvacancyb=.false.
test ${lvacancy} = 1 && lvacancyb=.true.
test ${lvacancy2} = 0 && lvacancy2b=.false.
test ${lvacancy2} = 1 && lvacancy2b=.true.
test ${lvacancy3} = 0 && lvacancy3b=.false.
test ${lvacancy3} = 1 && lvacancy3b=.true.
test ${lvacancy4} = 0 && lvacancy4b=.false.
test ${lvacancy4} = 1 && lvacancy4b=.true.
test ${lvacancy5} = 0 && lvacancy5b=.false.
test ${lvacancy5} = 1 && lvacancy5b=.true.
test ${limpurity} = 0 && limpurityb=.false.
test ${limpurity} = 1 && limpurityb=.true.
test ${limpurity2} = 0 && limpurity2b=.false.
test ${limpurity2} = 1 && limpurity2b=.true.
test ${lrec} = 0 && lrecb=.false.
test ${lrec} = 1 && lrecb=.true.
test ${lcheckcubic} = 0 && lcheckcubicb=.false.
test ${lcheckcubic} = 1 && lcheckcubicb=.true.



#label="${nx}${ny}${nz}_${nkx_r}${nky_r}${lkshifted}_${pbcz}_${phi_type}_ef${ef}V${V}tf${tf}T${T}_mu${lfixmu_k}${other}"
#label="${nx}${ny}${nz}_${nkx_r}${nky_r}${lkshifted}_${pbcz}_${phi_type}_ef${ef}V${V}tf${tf}_v${lvacancy}${lvacancy4}${lvacancy3}${lvacancy5}${lvacancy2}${other}"
#label="${nx}${ny}${nz}_${nkx_r}${nky_r}${lkshifted}_${pbcz}_${phi_type}_ef${ef}V${V}tf${tf}_v${lvacancy}${lvacancy4}${lvacancy3}${lvacancy5}${lvacancy2}_ec${ec0}${other}"
#label="${nx}${ny}${nz}_${nkx_r}${nky_r}${lkshifted}_${pbcz}_${phi_type}_ef${ef}V${V}tf${tf}_d${ldisorder_r}${ldisorder_r2}_c${conc_vac}_ec${ec0}_mu${lfixmu_kr}${lfixmu}${other}"
#label="9915_110_0_1-2_ef-6V3.5tf0.4_v10001"
#label="9915_110_1_1-2_ef-6V3.5tf0.4_v00100"
#label="9915_110_1_1-2_ef-0.1V3.5tf0.2_v00100"
#label="9915_110_0_1-2_ef-0.1V3.5tf0.2_v10000"
#label="9915_110_0_1-2_ef-0.1V3.5tf0.2_v10001"
#label="7715_110_1_1-2_ef-6V3.5tf0.4_v10001"
#label="7715_110_1_1-2_ef-0.1V3.5tf0.2_v10001"
#label="7715_110_1_1-2_ef-0.1V3.5tf0.2_v10010"
#label="5515_110_0_1-2_ef-0.5V3.5tf0.4_v10000"
#label="7715_110_1_1-2_ef-0.1V3.5tf0.2_v10010"
#label="7715_110_1_1-2_ef-0.1V3.5tf0.2_v10001"
#label="9915_110_1_1-2_ef-0.1V3.5tf0.2_v00100"
#label="7715_110_1_1-2_ef-6V3.5tf0.4_v10001"
#label="111115_110_0_1-2_ef-6V3.5tf0.4_v10000"
#label="test105"
#label="111115_110_0_1-2_ef-0.1V3.5tf0.2_v10000"
#label="111112_221_0_1-2_ef-0.5V3.5tf0.2_v10000"
#label="151510_110_0_kim_ef-6V3.3tf0.4_d10_c0.05_ec0_mu11"

cat > par.in << EOF
&input_par
loopk=${loopkb}
lscfk=${lscfkb}
lscfkr=${lscfkrb}
loopkr=${loopkrb}
loopr=${looprb}
lscfr=${lscfrb}
lscfloop=${lscfloopb}
lscfloop_kr=${lscfloop_krb}
lnscfr=${lnscfrb}
lpathnscfr=${lpathnscfrb}
lpathnscf=${lpathnscfb}
label = $label
ef_=${ef}
ec_=${ec}
tf_=${tf}
tc_=${tc}
b_init=${b_init}
Vcf_=$V
lambda_init=${lambda_init}
mu_init=${muinit}
ef0=$ef0
ec0=$ec0
nx=$nx
ny=$ny
nx_kr=${nx_kr}
ny_kr=${ny_kr}
nz=$nz
nkx_r=${nkx_r}
nky_r=${nky_r}
phi_type=${phi_type}
lvacancy=${lvacancyb}
lvacancy2=${lvacancy2b}
lvacancy3=${lvacancy3b}
lvacancy4=${lvacancy4b}
lvacancy5=${lvacancy5b}
limpurity=${limpurityb}
limpurity2=${limpurity2b}
ldisorder_r=${ldisorder_rb}
ldisorder_r2=${ldisorder_r2b}
conc_vac=${conc_vac}
T_min=$T
T_max=${T_max}
delta_T=${delta_t}
lreadparam=${lreadparamb}
lfixmu_kr=${lfixmu_krb}
lfixmu_k=${lfixmu_kb}
lfixmu=${lfixmu_b}
lmucf=${lmucfb}
ltr=${ltrb}
lkshifted=${lkshiftedb}
lkred_r=${lkred_rb}
pbcz=${pbczb}
T_g=${T_g}
eta_7=${eta_7} 
eta2_f7=${eta2_f7} 
eta3_f7=${eta3_f7}
eta1_d=${eta1_d}
eta2_d=${eta2_d}
eta3_d=${eta3_d}
lrec=${lrecb}
d_ec=${d_ec}
d_ef=${d_ef}
d_ec_kr=${d_ec_kr}
d_ef_kr=${d_ef_kr}
t_d1=${t_d1}
t_f1=${t_f1}
lfullbz=${lfullbzb}
d_e_rec=${d_e_rec}
alfa_nz=${alfa_nz}
eta_v=${eta_v}
eta_v2=${eta_v2}
lcheckcubic=${lcheckcubicb}
eta_v3=${eta_v3}
lambda_par_k=${lambda_par_k}
alpha_mu_k=${alpha_mu_k}
lk3d=${lk3d}
alpha_min_k=${alpha_min_k}
alpha_max_k=${alpha_max_k}
threshold=${threshold}
threshold_kr=${threshold_kr}
threshold_k=${threshold_k}
scf_type=${scf_type}
t_stop=${t_stop}
n_step=${n_step}
tolerance_rksuite=${tol_rk}
threshold_rksuite=${thr_rk}
time_integr=${time_integr}
period_time_ev=${period}
 /

&input_calc
nz_scf=${nz_scf}
lkshifted_r_nscf=${lkshifted_r_nscfb}
lwfc=${lwfc}
e_thr1=${e_thr1}
e_thr2=${e_thr2}
nkx_r_nscf=${nkx_r_nscf}
nky_r_nscf=${nky_r_nscf}
lbinning=${lbinningb}
nkpath_plot=${nkpath_plot}
lrhoktot=${lrhoktotb}
nz_plot=${nz_plot}
nkx_plot=${nkx_plot}
nky_plot=${nkx_plot}
lgreen=${lgreenb}
lborn=${lbornb}
lg0=${lg0b}
lv_nscf=${lv_nscfb}
omega_min=${omega_min}
omega_max=${omega_max}
delta_omega=${delta_omega}
size_h_kr_red_e=2000
en_min_qpi=-10
en_max_qpi=10
stm_layer_ratio=${stm_layer_ratio}
nkx_chern=${nkx_chern}
iplane_xyz=${iplane_xyz}
nk_dos=${nk_dos}
T_dos=${T_dos}
/
EOF

cp par.in  par_${label}.in

./tki_ga < par.in #| tee par_$label.out

#mv par_${label}.in ./$label/
#rm par_${label}.in
#echo  -ne ${d_ef}"\t"  >> $outfile
#awk  '{printf "%s\t", $1}' $label/dissipation >> $outfile
#grep "4   4   0"  $label/b_r2 | awk '{print $4, $5}' >> $outfile
#grep "4   4   0"  $label/b_r2 | awk '{print $4, $5}'

#awk 'FNR==1{ print $1}' $label/b_k2 > b_1
#awk 'FNR==1{ print $1}' $label/b_k2 
#export b_init="`cat b_1`"
#awk 'FNR==1{ print $2}' $label/b_k2 > lambda_1
#awk 'FNR==1{ print $2}' $label/b_k2 
#export lambda_init="`cat lambda_1`"
#awk 'FNR==1{ print $3}' $label/b_k2 > mu_1
#awk 'FNR==1{ print $3}' $label/b_k2 
#export muinit="`cat mu_1`"

#done
#done

done
