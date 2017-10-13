!include 'rksuite.f90'
include 'rkf45_new.f90'
include 'rksuite_90.f90'
!C4W95C
program tki
!  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R
use rksuite_90

IMPLICIT NONE

INTEGER, PARAMETER              :: idp=KIND(1.0D0), ids=KIND(1.0)
INTEGER, PARAMETER              :: idpc=kind((1.0d0,0.0d0)),idsc=kind((1.0,0.0))
REAL(KIND=idp), PARAMETER       :: Pi = 3.141592653589793D0
complex(KIND=idpc), allocatable	:: Ham(:,:),  W(:,:), V(:,:),Ham_k(:,:),  W_k(:,:),W_k_tot(:,:,:), W_k_tot_hsp(:,:,:),e_k_tot_hsp(:,:),&
				   V_k(:,:),Ham_k_red(:,:), w_k_red(:,:),rot_kdp(:,:),&
				   Ham_kr(:,:), W_kr(:,:),W_kr_rec(:,:),V_kr(:,:), W_kr_tot(:,:,:),w_kr_p(:,:),w_kr_m(:,:), ham_kr_p(:,:), ham_kr_m(:,:), ham_kr_pm(:,:),ham_kr_mp(:,:),&
				   phi_sigma_alpha(:,:,:,:,:), phi_alpha_sigma(:,:,:,:,:),V_k_tot(:,:,:),v_kr_tot(:,:,:),kin_energy(:,:,:,:,:,:),kin_energy_k_tot(:,:,:,:,:,:),&
				   kin_energy2(:,:,:,:,:,:),kin_energy2_k_tot(:,:,:,:,:,:),kin_energy2_k_nscf(:,:,:,:,:,:),kin_energy2_k_path3d(:,:,:,:,:,:),&
				   H_psi(:),H_psi_tot(:,:,:),H_kr_psi_tot(:,:,:),kin_energy_k_nscf(:,:,:,:,:,:),kin_energy_k_path3d(:,:,:,:,:,:),&
				   phi_k_sigma_alpha_tot(:,:,:,:,:),phi_k_alpha_sigma_tot(:,:,:,:,:),phi_kr_sigma_alpha_tot(:,:,:,:,:),phi_kr_alpha_sigma_tot(:,:,:,:,:),kin_energy_kr_tot(:,:,:,:,:,:),&
				   phi_k_sigma_alpha_nscf(:,:,:,:,:),phi_k_alpha_sigma_nscf(:,:,:,:,:),phi_kr_sigma_alpha_nscf(:,:,:,:,:),phi_kr_alpha_sigma_nscf(:,:,:,:,:),&
				   corr_fct_kr_nscf(:,:), corr_fct_k_nscf(:,:), V_kr_nscf(:,:), V_k_nscf(:,:),w_kr_tot_nscf(:,:,:), w_r_tot(:,:,:),kin_energy_kr_nscf(:,:,:,:,:,:),&
   				   phi_k_sigma_alpha_path3d(:,:,:,:,:),phi_k_alpha_sigma_path3d(:,:,:,:,:),H_kr_tot_nscf(:,:,:),kin_energy2_kr_tot(:,:,:,:,:,:,:),kin_energy2_kr_nscf(:,:,:,:,:,:,:),&
   				   kin_energy3_kr_tot(:,:,:,:,:,:,:),kin_energy3_kr_nscf(:,:,:,:,:,:,:),kin_energy3(:,:,:,:,:,:),&
   				   phi_sigma_alpha2(:,:,:,:,:), phi_alpha_sigma2(:,:,:,:,:),phi_k_sigma_alpha_hsp(:,:,:,:,:),phi_k_alpha_sigma_hsp(:,:,:,:,:),&
   				   phi_kr_sigma_alpha2_tot(:,:,:,:,:,:),phi_kr_alpha_sigma2_tot(:,:,:,:,:,:),phi_kr_sigma_alpha2_nscf(:,:,:,:,:,:),phi_kr_alpha_sigma2_nscf(:,:,:,:,:,:),&
   				   phi_sigma_alpha3(:,:,:,:,:), phi_alpha_sigma3(:,:,:,:,:),kin_energy_k_hsp(:,:,:,:,:,:),&
   				   phi_kr_sigma_alpha3_tot(:,:,:,:,:,:),phi_kr_alpha_sigma3_tot(:,:,:,:,:,:),phi_kr_sigma_alpha3_nscf(:,:,:,:,:,:),phi_kr_alpha_sigma3_nscf(:,:,:,:,:,:),&
   				   phi_k_sigma_alpha_chern(:,:,:,:,:),phi_k_alpha_sigma_chern(:,:,:,:,:),kin_energy_k_chern(:,:,:,:,:,:),e_kr_cmplx(:),&
					h_psi_tot0_k(:,:,:),h_psi_tot1_k(:,:,:)
complex(KIND=idpc), allocatable	:: U(:,:),Um1(:,:), Ham_rot(:,:), UstarU(:,:),Ham_input(:,:,:,:,:), Ug78(:,:), Ug78m1(:,:),Ham_input_temp(:,:,:,:,:), symm_hsp(:,:),eigen_symm(:),w_symm(:,:),&
				   w_symm_xy_pe(:,:),w_symm_xy_pf(:,:),w_symm_xy_me(:,:),w_symm_xy_mf(:,:), symm_hsp_p(:,:),symm_hsp_m(:,:),symm_hsp_pm(:,:),symm_hsp_mp(:,:), simm_xy(:,:),M_xy(:,:),&
				   M_xy_tmp_f(:,:),M_xy_tmp_e(:,:),M_xy_tmp_f2(:,:), m_check(:,:), m_check2(:,:), m_xzyz(:,:), m_x(:,:), m_y(:,:),M_yz(:,:),w_xz(:,:),&
				   w_x_p(:,:),w_x_m(:,:),w_y_p(:,:),w_y_m(:,:),m_x_kr(:,:), w_m_x_kr(:,:),m_x_kr_pi(:,:),w_m_x_kr_pi(:,:),w_x_m_kr_pi(:,:),w_x_p_kr_pi(:,:),&
				   w_x_p_kr(:,:),w_x_m_kr(:,:),w_y_p_kr(:,:),w_y_m_kr(:,:),w_xy_p_kr(:,:),w_xy_m_kr(:,:),W_X_KR_PI(:,:),W_X_KR(:,:),w_yz_m_kr(:,:),w_yz_p_kr(:,:),&
				   HAM_KR_TMP_P(:,:),HAM_KR_TMP_m(:,:),w_x(:,:), w_y(:,:),HAM_KR_TMP_pm(:,:),HAM_KR_TMP_mp(:,:),w_kr_p_tmp(:,:), w_kr_m_tmp(:,:), u_symm_kr_p(:,:), u_symm_kr_m(:,:),&
				   e_kr_tot_nscf_cmplx(:,:),w_k_tot_nscf(:,:,:),Ham1k(:,:,:),Ham0k(:,:,:)

				   
real(KIND=idp), allocatable	:: E(:),E_k(:),E_k_tot(:,:),E_kr(:),E_kr_tot(:,:), nf(:),nc(:), ec(:), ef(:,:),  nfc(:), b(:),  b_new(:),&
				   lambda(:),  lambda_new(:),nf_kr(:),nff_kr(:),nfk_kr(:), nc_kr(:), ec_kr(:), ef_kr(:), nfc_kr(:),nfck_kr(:), b_kr(:), b_new_kr(:), lambda_kr(:), lambda_new_kr(:), &
				   one(:), zero(:),one_kr(:), zero_kr(:), alpha_jzi(:,:,:),UtU(:,:),alpha_jzi_mat(:,:), alpha_msigma(:,:,:,:), dos(:,:),  nff(:), &
				   nn_theta_phi(:,:), nn_xyz(:,:),nn2_theta_phi(:,:), nn2_xyz(:,:),nn3_theta_phi(:,:), nn3_xyz(:,:),nff_imag(:),nfc_imag(:),&
				    k_vecs_nscf(:,:),delta_k_nscf(:),delta_k_path3d(:),delta_k_hsp(:),delta_k_chern(:),e_k_tot_nscf(:,:),ew_k_tot_nscf(:,:),delta_kr_tot(:),E_r_tot(:,:),&
				   kr_vecs_nscf(:,:),delta_kr_nscf(:),e_kr_tot_nscf(:,:),ew_kr_tot_nscf(:,:,:,:),fermien_k(:,:),occup_k(:),occup_kr(:),fermien_kr(:,:),kr_vecs_r(:,:),&
				   e_k_tot_path3d(:,:),ew_k_tot_path3d(:,:,:),k_vecs_path3d(:,:),k_vecs_hsp(:,:),kr_vecs_nscf_klarge(:,:),kr_vecs_nscf_qpi(:,:),&
				   e_kr_tot_nscf_pm(:,:,:),w_kr_tot_nscf_pm(:,:,:,:),e_kr_p(:), e_kr_m(:),e_x(:), e_y(:),e_m_x_kr(:),e_m_x_kr_pi(:),e_x_kr(:),e_x_kr_pi(:),&
				   ew_kr_tot_p_nscf(:,:,:,:),ew_kr_tot_m_nscf(:,:,:,:),spin_kr_tot_nscf(:,:,:,:,:),spin_kr_tot_m_nscf(:,:,:,:,:),spin_kr_tot_p_nscf(:,:,:,:,:),&
				   ew_kr_tot_nscf_rec(:,:,:,:),spin_kr_tot_nscf_rot(:,:,:,:,:),spin_kr_tot_p_nscf_rot(:,:,:,:,:),spin_kr_tot_m_nscf_rot(:,:,:,:,:),e_k_tot_path3d_red(:,:), e_k_red(:),&
				   !k_vecs(:,:),
				   k_vecs_chern(:,:),parity_hsp_f(:,:),parity_hsp_e(:,:),&
				   k_vecs_kr(:,:),dos_kr(:,:),k_vecs3d(:,:),k_vecs2d(:,:), k_vecs_r(:,:), factor_k_tot(:),factor_kr_tot(:), delta_k_tot(:),factor_k_nscf(:),&
				   factor_kr_nscf(:),factor_k_path3d(:),eigen_symm_real(:),&
				   ec_site(:,:), ef_site(:,:),ef_site_k(:),ec_site_k(:), e_corr_fct_kr_nscf(:,:),e_corr_fct_k_nscf(:,:), &
				   k_vecs_r_rot(:,:),k_tot(:,:),ef_site_kr(:,:),ec_site_kr(:,:),k_vecs_r_arpes(:,:),ef_site0(:,:),ec_site0(:,:),ef_site1(:,:),ec_site1(:,:),ef0_site(:,:)
integer				:: nkx, nky, nkz, Nktot, ixyz,ixyz2, k_total_size, n_k, size_h_kr, n_k_kr, n_sites_kr, nkx_plot,nky_plot,nkz_plot,nk_plot,n_kvecs_irrid, nk3d_irr, nk2d_irr,nkr_plot, &
				   nkpath_plot, nkpathr_plot, nkx_r, nky_r,nk_plot_path3d,nk_plot_hsp, nk_r,nk_r_nscf,nkx_r_nscf,nky_r_nscf,nkpath_r, nk_r_rot, n_k_total,nkpath_x,nkpath_y,&
				   nkpath_xy, nk_x_rot,  nk_y_rot, n_shells, n_shells_in,size_h_kr_red_e,nkxy_plot,nkx_kr, nky_kr,&
				   nkpath_plot_xy, nk_plane,nkpath_r_xy, nkx_kk, nky_kk , nkx_k, nky_k, nk_r_tot, nk_r_arpes, ik_arpes, nx_g, ny_g, nkx_l, nky_l, ysameasx, sx, sy, swap_xy, nz_scf, nz_g0,&
				   nk_chern, nkx_chern,nky_chern,iplane, nkxy_chern, nkz_chern
integer, allocatable		:: nn_ijk(:,:),nn2_ijk(:,:),nn3_ijk(:,:), index_psi(:,:),index_psi_tot(:,:,:), index_kr_psi_tot(:,:,:), list_kk(:,:),nn_ss(:,:),H_t(:)
integer				:: nx, ny, nz, Ntot,total_size, n_sites,  nk3d, nk2d,iz_ind,size_h_kr_red, l_orbital,l_orbital_f, l_orb(0:1), l_orb_tot, nx_kr, ny_kr
integer				:: ix,iy,iz,icf,is, o1, o2,o3,k2,is3,is4, dim_hk,n,lx,ly,lz,nkr_plot_klarge,nkxpath_plot,nkypath_plot,n_sym_arpes,nkxypath_plot	!x,y,z, c or f, up /down +/- 
integer				:: ix2,iy2,iz2,icf2,icf3,is2	!x,y,z, c or f, up /down +/- 
integer				:: ikx,iky,ikz, ik!,icf,is	!x,y,z, c or f, up /down +/- 
integer				:: ikx2,iky2,ikz2, ik2,ik_px,ik_py!,icf2,is2	!x,y,z, c or f, up /down +/- 
integer				:: N_nn, n_nn2,n_nn3		!nearest neighbours
integer				:: i,j,k,l,m, n_occ, n_empty,ikn
integer		:: ind, ind2, indp, ind2p, ind3, ind3p,ind4, ind4p,ind5p, ind_vac, ind_vac2, ind_vac3, ind_vac4, ind_vac5, ind1,ind5,ind6,ind7,ind8,ind9,ind10
integer		::i_iter, max_iter, i_iter_mu, n_ksf, nx_red, ny_red
logical 	::ok, pbcx,pbcy,pbcz, Vtki, lkx,lky,lkz,lkx2,lky2,lkz2,lkx3,lky3,lkz3, &
		  loopk,loopr, loopkr,ldos, lfixmu,lnscfk,lwfc, lysameasx,lzsameasx, lkred, lpathnscf, lscfk, lscfkr, lreadbk, lreadbkr, lfileexists,&
		  lreadparam, lb2tf, lreadbr, lcorr_fct, lrednscf, lspectrfunct, lscfr, lnscfr, lreadwfc, lvacancy, limpurity, limpurity2,lpathnscfr, lwritewfc,lvacancy_kr,lvacancy_kr2,lvacancy_kr3,&
		  lsfrkred,lvacancy2,lvacancy3, lfullbz, ldisorder_r,ldisorder_r2, ldisorder_r_bulk, lvacancy_kr4, lvacancy_kr5, lvacancy_kr6,limpurity_kr,limpurity_kr2, lkshifted_r,&
		  lkshifted, ltr,lvacancy4,lvacancy5, limpurity_z, limpurity_x, limpurity_i,lfft,lqpi,ljdos,larpes,lphaseshifts,lmagimp,l_k_dz,l7z,&
		  lspectrfunct_kr, lkred_r, lfixmu_kr,lfixmu_k, lbinning, lwfctot, lrhoktot ,lkshifted_r_nscf, lgreen,lborn, lg0, lscfloop, lscfloop_kr, lv_nscf, lvkn, lmucf, lhgte, lram,ldiag_g0k,&
		  lreadham, l2dmodel, lphioverlap, lkinoverlap, lrec, lc4v, lc2v, lcs,lnos, lfftreal, lmuscf,lmufscf,lcheckcubic,lchern
real(KIND=idp)	::theta,phi, v_1, v_2, omega_cutoff, theta_xy,lambda_dec_kr, norm,norm2 ,det_rotmat, r2_temp(3), t_dos, nf_0,nc_0,t_stop,period_time_ev
integer		:: r1a,r2a,r3a, r1b,r2b,r3b, r1c,r2c,r3c, det_r1r2r3,r1r2r3(3,3), r1r2r3_inv(3,3), planeindex(3),r1dotr2, nk_dos,n_step
real(KIND=idp)	::e_max, e_min,delta_e, energy, T_min, T_max, delta_T, delta_k, delta_kr, gamma_lor, thr_fermi, ef0,ec0, occup_tot, delta_ef,N_el, e_thr1, e_thr2, T_g, ef0_up, ef0_down, ec0_up, ec0_down,&
		  d_ec, d_ef7, d_ef, d_ec_kr, d_ef7_kr, d_ef_kr, d_ef_kr_82, en_min_qpi, en_max_qpi, d_e_rec, alfa_nz, d_ec_2_kr, d_ef7_2_kr, d_ef_2_kr,d_e_rec2,k1_c,k2_c,k3_c,delta_c,delta_time
integer		:: nxn, nyn,neqn
 character(LEN=100)::label
 character(LEN=3)::phi_type, surf_type
 character(LEN=2)::orb_type,scf_type
real(KIND=idp)			::h1,f1,h2,f2,hx,h7,f7,m78,kz,kz2,kpar2,h72,g1d,g2d,l1d,l2d,g1f,g2f,l1f,l2f,g7,l7,k2xy, eckp1,eckp2,efkp1,efkp2,efkp7,td12,tf12,b_kr0,b_kr1,en,total_weight,f1k2,h1k2
real(KIND=idp)		:: ec_, ef_, tc_, tf_, Vcf_, filling, b_init, lambda_init, mu_init, eta_v, ef_7, eta_vz, etav2_x1,etav2_x2, etav2_z1,etav2_z2,etav2_x7,etav2_z7,&
				 etav3_x1,etav3_x2, etav3_z1,etav3_z2,etav3_x7,etav3_z7, etav3_x1b,etav3_x2b, etav3_z2b, real_h, imag_h, ef_2, ec_2, stm_layer_ratio,weight
real(KIND=idp)		:: ec_n, ef_n, tf_n, Vcf_n, T_K, eta, eta_d, eta_f, eta_df(0:1), eta2_d, eta2_f7,eta2_f, eta2_f78,eta2_df(0:1),eta2_ok, eta1_d, eta1_f, eta1_df(0:1), &
			   eta_7, eta_v7, eta_78, eta3_f7, eta3_f8, eta3_78, eta3_ok, eta2_dx, eta2_f81,eta2_f82,eta2_f812, eta2_f781, eta2_f782, eta3_d, eta_v2, eta_v3
real(kind=idp)		:: T, mu,  mu_new,mu_old, alpha,alpha_k,alpha_kr,alpha_k_u,alpha_min_k, alpha_max_k, alpha_mu, threshold,threshold_k,threshold_kr, &
			   diff, diff_mu, thr_mu, mu_guess, dndmu, lambda_par,lambda_par_k,&
			   lambda_par_kr,alpha_min, alpha_max, delta_alpha, alpha_min_lambda, alpha_max_lambda, alpha_lambda, diff_thr
real(kind=idp)		:: Nftot, Nctot, Nfctot, fctot, E_tot, E_tot_l, E_tot_k, E_tot_k_l, N_el_k_c
real(kind=idp)		:: b_k, lambda_k, b_new_k, lambda_new_k, mu_k, mu_new_k,mu_new_kr, f_k,f_kr, mu_old_k, alpha_mu_k,alpha_mu_k_u,alpha_mu_kr,nf_k,&
			   nfc_k,nc_k,ec_k,ef_k,nff_k, mu_kr, occup,nf7_k, occupation_k, conc_vac, delta_e7_read, delta_e8_read, delta_tf_read, delta_v_read,delta_tc_read
real(kind=idp)			:: n_el_k,n_el_kr,threshold_rksuite, tolerance_rksuite
integer			:: numtd
real(KIND=idp)	::time_0, time_1, time_init, time_k, time_buildh, time_diag, time_nf, time_mu, time_sc, time_wr,time_buildh_k,time_diag_k,&
		  time_wrk, time_nfk, time_setvark, time_zhpev, time0_zhpev,time1_zhpev, time_kr, time_r, time_kr_setvar, time_kr_buildh, time_kr_diag, time_kr_nf, time_kr_mu, time_kr_write,&
		  time_nscfk, time_nscfkr, time_sfr, time_rnscf, time_invert, time_0_inv, time_1_inv,time_matmul, time_0_matmul, time_1_matmul, time_write, time_fourier, time_vkn, time_rhovrho,&
		  time_green,time_green_0,time_green_1,time_r8f2,time_time_ev,time_r8f2_0,time_time_ev_0,time_nfimag,time_nfimag_0,time_E, time_e_0
integer		::n_max_link, nz_plot,n_states_in,n_entry_in,j1, n_link_tot, iplane_xyz
integer, allocatable	::z_values_plot(:),list_k_vecs_r_arpes(:,:)
 character(len=80):: format_smear, format_wr, format_er, format_k
 logical, allocatable	:: lvac_r(:), lvac_kr(:)
complex(kind=idpc):: weight_cmplx
integer		:: n_vac_kr, n_imp_kr, n_vac_r, n_imp_r , dist, lk3d
real(kind=idp):: omega_min, omega_max, delta_omega, omega, t_d1, t_d2, t_f1, t_f2, t_f3, dist_max,  j3
!integer, allocatable		:: z_vac_kr(:), z_imp_kr(:),x_vac_r(:,:), x_imp_r(:,:)
complex(KIND=idpc),allocatable::sigma_xyz_cf(:,:,:,:,:,:), sigma_pauli(:,:,:), simm_refl(:,:,:,:,:,:),simm_spin(:,:,:),simm_orb(:,:,:),simm_orb_d(:,:,:), simm_oper(:,:,:),simm_oper2(:,:,:)
 character(len=80):: format_in, format_out, label_in, label_out, format_p, format_h


complex(KIND=idpc),allocatable:: w_k_tot_chern_fp(:, :,:),w_k_tot_chern_fm(:,:,:),w_k_tot_chern_ep(:, :,:),w_k_tot_chern_em(:,:,:),e_k_tot_chern_p(:,:),e_k_tot_chern_m(:,:),&
				 M_xy_f(:,:),M_xy_e(:,:),w_k_f(:,:),w_k_e(:,:),w_xy_f_tmp(:,:),w_xy_e_tmp(:,:),w_xy_f(:,:),w_xy_e(:,:),&
				 u_k_chern(:,:,:),u_mat_fp(:,:),u_mat_fm(:,:),u_mat_ep(:,:),u_mat_em(:,:), w_u_em(:,:),w_u_ep(:,:), w_u_fm(:,:), w_u_fp(:,:),&
				 e_u_em(:),e_u_ep(:), e_u_fm(:), e_u_fp(:), f12_k(:,:), det_(:), ham_k_p(:,:), ham_k_m(:,:),ham_k_pm(:,:),ham_k_mp(:,:),w_k_p(:,:), w_k_m(:,:),&
				 w_xy(:,:),w_xy_p(:,:),w_xy_m(:,:), ham_k_tmp(:,:), w_k_p_tmp(:,:), w_k_m_tmp(:,:),w_xzyz_p(:,:),w_xzyz_m(:,:),w_xzyz(:,:), w_m_p(:,:), w_m_m(:,:),&
				 w_yz_m(:,:),w_yz_p(:,:),w_yz(:,:),g_k(:,:), g_k_res(:,:),g_k2(:,:)
real(KIND=idp),allocatable:: e_xy_f(:),e_xy_e(:), e_k_p(:),e_k_m(:),e_xy(:),e_xzyz(:),e_yz(:), dos_k(:)
real(KIND=idp):: alfa_rot_mat(3,3), beta_rot_mat(3,3), gamma_rot_mat(3,3), alfa_rot, beta_rot, gamma_rot, rot_mat_surf(3,3), rot_mat_temp(3,3), k_surf(2,2), x_surf(2,2),rot_mat_surf_check(3,3),&
		r1square,r2square,cosr1r2
real(kind=idp),allocatable:: b_kr_0(:),lambda_kr_0(:)
 character(len=3):: time_integr



numtd=4
call omp_set_num_threads(numtd)



NAMELIST /input_par/ 	nx,ny,nz, label, b_init, lambda_init, mu_init, phi_type,lfixmu, ef_,ec_, Vcf_, tf_, ef0, ec0, lreadparam, lkshifted,ef0_up, ec0_up, ef0_down, ec0_down, &
			pbcx,pbcy,pbcz, nkx_r, nky_r, nkx_kk, nky_kk,lkred,lkred_r,T_min,lfixmu_kr,lfixmu_k,T_max, delta_T,l_k_dz,ef_7, eta_7, eta_v7, eta_78,&
			lvacancy,lvacancy2,lvacancy3,lvacancy_kr, lvacancy_kr2, lvacancy_kr3, lvacancy_kr4, lvacancy_kr5,ldisorder_r,ldisorder_r2,conc_vac,limpurity_kr,limpurity_kr2,&
			lvacancy4,lvacancy5, ltr, limpurity, limpurity2, limpurity_z, limpurity_x, limpurity_i,lmagimp,ldisorder_r_bulk,diff_thr,&
			loopk, loopkr, loopr,lscfr, lnscfr,lpathnscfr,lpathnscf, t_g, lscfloop, lscfloop_kr,lmucf, lscfk, eta,eta_f, lscfkr, tc_, eta_v,lambda_par_k, eta2_d, eta2_f, eta1_d, eta1_f,&
			eta2_f7, eta2_f78, eta3_f7, eta3_f8, eta3_78,eta_vz,eta2_dx,eta2_f81,eta2_f82,eta2_f812, eta2_f781, eta2_f782, eta3_d,&
			etav2_x1,etav2_x2, etav2_z1,etav2_z2,etav2_x7,etav2_z7,eta_v2, eta_v3,&
			etav3_x1,etav3_x2, etav3_z2,etav3_x1b,etav3_x2b,etav3_z2b,etav3_x7,etav3_z7,nf_0,&
			t_d1, t_d2, t_f1, t_f2, t_f3, lreadham, dist, n_shells,lfullbz,alpha_mu_k,lmuscf,lmufscf,surf_type,lk3d,&
			d_ec, d_ef7, d_ef, d_ec_kr, d_ef7_kr, d_ef_kr,d_ef_kr_82, delta_e7_read, delta_e8_read, delta_tf_read,delta_v_read,delta_tc_read,lcheckcubic,&
			v_1,v_2, omega_cutoff, theta_xy, ef_2, ec_2, nx_kr, ny_kr, lrec, eta_d, eta2_d, eta3_d, d_e_rec, d_e_rec2,alfa_nz, d_ec_2_kr, d_ef7_2_kr, d_ef_2_kr,lambda_dec_kr,&
			alpha_min_k, alpha_max_k,lchern,threshold,threshold_k,threshold_kr,scf_type,t_stop,n_step,threshold_rksuite, tolerance_rksuite,time_integr,period_time_ev
NAMELIST /input_calc/ 	nkpath_plot, nkx_plot, nky_plot, nkx_r_nscf, nky_r_nscf,lspectrfunct, lwfc, e_thr1, e_thr2,lspectrfunct_kr, lbinning, lrhoktot,lkshifted_r_nscf, nz_plot, nx_g, ny_g, lgreen,&
			omega_min, omega_max, delta_omega,lborn,lv_nscf, lvkn, nz_scf, larpes,lphaseshifts, lg0, en_min_qpi, en_max_qpi,size_h_kr_red_e,stm_layer_ratio, ljdos, nkx_chern, iplane_xyz,&
			r1a,r1b,r1c,r2a,r2b,r2c,r3a,r3b,r3c,b_kr0,b_kr1,nk_dos,t_dos

b_kr0=1
b_kr1=1
nf_0=0.5
lmufscf=.false.
lchern=.false.
l7z=0
scf_type='sb'
threshold_rksuite=1e-7
tolerance_rksuite=1e-7
time_integr='rks'


surf_type="001"

orb_type="f5"

n_shells_in=3
n_states_in=18
format_in="(5i5,2f16.10)"
n_entry_in=(2*n_shells_in+1)**3*n_states_in**2

!label="8830V314_f02_ppo_e-55_phisz"
!label="5530_12_-01_nomu"
!label="test"
!nx=3
!ny=3
!nz=20
nkx_kk=1	!number of periodic replicas in k-loop (to better fix mu)
nky_kk=1
!nkx_r=1	!number of periodic replicas in r-loop
!nky_r=1
!
!pbcz=0
!phi_type="1-2"		!1-2, 3-2, 5-2, g_1, g_2, g_3, dze, ide, six, trv
!
nkpath_plot=50
!nkpathr_plot=10
nkx_plot=50
nky_plot=50
!nkx_r_nscf=6
!nky_r_nscf=6
lspectrfunct=1
lspectrfunct_kr=1
nkx_chern=101
nky_chern=101

!
!lkshifted=0	!in k and kr loop, only along xy
lkred=1
ltr=1		!time reversal for generation k points in r loop (not reduced routine)
lysameasx=1	! if lkred=1
lzsameasx=0	! if lkred=1
!
!
lmucf=1
lmuscf=1	!sets mu to the gap
lreadparam=1	!size, energies ef_, hoppings
loopk=1	!k loop
 lscfk=1	!k scf, and writes to file b_k2
 lreadbk=1	!reads b_k from file b_k2
loopkr=1	!kr loop
 lscfkr=1	!kr scf, and writes to file b_kr2
 lreadbkr=1	!reads b_kr from file b_kr2
 lnscfk=1	!does a nscf run
 lpathnscf=1	!does a nscf run along a path (if not a nscf run)
loopr=1
 lreadbr=1	!reads b lambda mu , ef ec
 lscfr=1
 lnscfr=1
 lpathnscfr=1
 lsfrkred=1	!use the reduced k-K procedure for lpath (only when nx=ny)			=0 for bands in reduced zone
 lfullbz=1	!for lpathnscfr=1 generates k points in the full bz, or in the reduced one	=0 for bands in reduced zone and now for kr path too lpathnscf
 lreadwfc=0	!reads wfc
 lwritewfc=0	!writes wfc
 lrhoktot=0	!1=keeps al k in memory for NSCF; 0=writes each k on disk, but not sum (has to be done by smear_r_nxnynz_k.sh)
 lcheckcubic=1
!
lvacancy=0	!vacancy in middle of top layer
lvacancy2=0	!vacancy in middle of bottom layer
lvacancy3=0	!vacancy in middle of middle layer
lvacancy4=0	!vacancy in middle of second layer
lvacancy5=0	!vacancy in middle of lastsecond layer
ldisorder_r=0
ldisorder_r2=0
ldisorder_r_bulk=0
conc_vac=0.5
lvacancy_kr=0	!vacancy on first layer kr
lvacancy_kr2=0	!vacancy on last layer kr
lvacancy_kr3=0	!vacancy in the middle layer kr
lvacancy_kr4=0	!vacancy on second layer kr
lvacancy_kr5=0 !vacancy in the last second layer kr
limpurity_kr=1	!impurity on first layer kr
limpurity_kr2=1	!impurity on last layer kr
!
limpurity=0
limpurity2=0
limpurity_z=0	!prop to tau_z
limpurity_x=0	!prop to tau_x
limpurity_i=0	!prop to tau_0=id
lmagimp=0
ec0=0
ef0=0
ec0_up=0
ef0_up=0
ec0_down=0
ef0_down=0
d_ec=0
d_ef=0
d_ef7=0
d_ec_kr=0
d_ef_kr=0
d_ef7_kr=0
d_ec_2_kr=0
d_ef7_2_kr=0
d_ef_2_kr=0
eta=0
eta_v=1
eta_vz=0
l_k_dz=.true.
eta2_d=0
eta2_f=0
eta2_f7=0
eta2_f78=0
eta1_d=0
eta1_f=0
!
l2dmodel=0
size_h_kr_red_e=2000
!

n_shells=3
dist=27
lreadham=0
t_d1=1
t_d2=1
t_f1=-0.1
t_f2=-0.1
t_f3=-0.1
en_min_qpi=-0.01
en_max_qpi=+0.01
!
delta_tf_read=1
delta_tc_read=1
delta_v_read=1
!
stm_layer_ratio=0.2
alfa_nz=0.2
lambda_dec_kr=0.0001
!
!ec_=0
!ef_=-0.1
!tf_=0.3
!Vcf_=4
filling=1.0000
!r
threshold=1e-2		!for the rSCF loop
lambda_par=0.12		!only alfa_lambda*lambda_par matters2
alpha_min=0.1		!mu and b, r
alpha_mu=0.1		!for the mu-loop, r
thr_mu=1e-8		!for computing mu (r loop)
!k
threshold_k=1e-7	!for the kSCF loop
lambda_par_k=0.1 !0.12			!lambda_par	!only alfa_lambda*lambda_par matters
alpha_min_k=0.4	!0.4			!alpha_min	!k loop
alpha_max_k=0.4 !0.4		!mu and b, r
alpha_mu_k=0.2			!alpha_min	!for the scf loop, there is no mu loop
diff_thr=10
!kr
threshold_kr=1e-3	!1e-3	!for the krSCF loop
lambda_par_kr=lambda_par_k!*0.1
alpha_kr=alpha_min!*0.1	!kr loop
alpha_mu_kr=0.01		!for the mu-loop, or for mu
!
!
T_min=0.01d0
T_max=5.500001d0
delta_T=0.01d0
T_g=0.01d0	!broadening for green's functions
!
! not to change
tc_=-1
pbcx=1
pbcy=1
lb2tf=1		
lfixmu=1	!1=does not sc-fix mu, for r (keeps initial mu)
lfixmu_kr=1	!1=does not sc-fix mu, for kr(keeps initial mu)
lfixmu_k=0	!1=does not sc-fix mu, for k (keeps initial mu)
Vtki=1
lhgte=0
!
nx_kr=1
ny_kr=1
!
lwfc=0
e_thr1=-0.1
e_thr2=0.1
lbinning=0
!
iplane_xyz=3 !compute chern number on k_{x,y,z}=0,pi plane

b_init=0.5
lambda_init=ef_-0.1

nz_scf=1		!must have nz_plot<=nz_scf

!ef_2=-532
!ec_2=-532

!!!!!!
READ(5,input_par)
size_h_kr=2*nz*(l_orbital+l_orbital_f)
size_h_kr_red_e=size_h_kr
nk_dos=2*(nz/2)
t_dos=t_min
READ(5,input_calc)
!!!!!!
!nz_plot=nz
!nkx_plot=nx

nz_g0=max(nz_scf, nz_plot)

nkxy_plot=nkx_plot
nky_plot=nkxy_plot/ny_kr
nkx_plot=nkxy_plot/nx_kr
if (nx_kr>1) nkx_plot=nkx_plot+1
if (surf_type=="110") nky_plot=int(nkx_plot/sqrt(2.0))
if (surf_type=="210") nky_plot=int(nkx_plot/sqrt(5.0))

nky_chern=nkx_chern
ef_2=ef_
ec_2=ec_
!if (abs(ef_2-532)<1e-5) ef_2=ef_
!if (abs(ec_2-532)<1e-5) ec_2=ec_


if (phi_type=="p33" .or. phi_type=="p31".or. phi_type=="pd4") then	!l=1 J=3/2
orb_type="p3"
j1=1
j3=1.5d0
elseif (phi_type=="p11".or.phi_type=="pd2") then !l=1 J=1/2
orb_type="p1"
j1=1
j3=0.5d0
else !l=3 J=5/2
orb_type="f5"
j1=3
j3=2.5d0
endif

l_orbital=1
l_orbital_f=1

if(phi_type=="q2d") then	!quasi 2d tb model non topological
lkinoverlap=0
lphioverlap=0
l_orbital=0
l_orbital_f=1
b_init=1
mu_init=0
lambda_init=0
lscfkr=0
lscfk=0
endif
if(phi_type=="sm2"  .or. phi_type=="trv" .or. phi_type=="ons") then
lkinoverlap=0
endif
if(phi_type=="smb") then
l_orbital=2 
l_orbital_f=2
endif
if(phi_type=="sm7") then
l_orbital=2
l_orbital_f=3
endif
if(phi_type=="dc2") then
l_orbital=0
l_orbital_f=1
nz=1
nx=1
ny=1
loopk=0
lscfkr=0
lscfloop_kr=0
lscfloop=0
lnscfr=0
lscfr=1
lscfloop=0
lscfk=0
b_init=1
mu_init=0
lambda_init=0
l2dmodel=1
endif

if(phi_type=="sra") then
tc_=-1
tf_=-tc_
ec_=-ef_
Vcf_=tc_	!<0
lscfloop=0
lscfloop_kr=0
lscfk=0
lscfkr=0
b_init=1
mu_init=0
lambda_init=0
l_orbital=1
l_orbital_f=1
endif

if(phi_type=="sxy") then
tc_=-1
!tf_=-tc_
!ec_=-ef_
lscfloop=0
lscfloop_kr=0
lscfk=0
lscfkr=0
b_init=1
mu_init=0
lambda_init=0
l_orbital=1
l_orbital_f=1
endif

if(lg0) ljdos=1

lphioverlap=0
lkinoverlap=0
if (phi_type=="1-2" .or. phi_type=="3-2" .or. phi_type=="5-2" .or. phi_type=="g_1" .or. phi_type=="g_2" .or. phi_type=="g_3" .or. phi_type=="g7_" .or. orb_type=="p3" .or. orb_type=="p1") then
!if (phi_type=="3-2" .or. phi_type=="5-2" .or. phi_type=="g_1" .or. phi_type=="g_2" .or. phi_type=="g_3" .or. phi_type=="g7_" .or. orb_type=="p3" .or. orb_type=="p1") then
lphioverlap=1
lkinoverlap=1
endif
if (phi_type=="1-2") lkinoverlap=0

l_orb(0)=l_orbital
l_orb(1)=l_orbital_f
l_orb_tot=l_orbital+l_orbital_f

dim_hk=2*(l_orbital+l_orbital_f)

if (loopr .and. lscfr .and. (lborn .or. lgreen .or. lg0)) lpathnscf=0
if (lpathnscf) nz_g0=nz
print *, "nz_g0=",nz_g0



eta_d=eta
!eta_f=eta
eta_df(0)=eta_d
eta_df(1)=eta_f
eta2_df(0)=eta2_d
eta2_df(1)=eta2_f
eta1_df(0)=eta1_d
eta1_df(1)=eta1_f

!print *, "Eta, Eta_1, Eta_2:"
!print *, eta_df
!print *, eta1_df
!print *, eta2_df

lmucf=1
if (t_max-T_min>0.01)  then
 lnscfk=0
 loopkr=0
 loopr=0
endif
!nkz_plot=nz
!lambda_init=-2.93
!mu_init=0
nkpath_r=nkpath_plot/nx
nkpath_plot_xy=int(sqrt(2.0D0)*nkpath_plot)
nkpath_r_xy=int(sqrt(2.0D0)*nkpath_r)
lkshifted_r=lkshifted	!only in r loop

T_K=-2*tc_*exp(-3.14*ef_*tc_/Vcf_**2)	!approximated tk

if(phi_type=="kyz") lysameasx=0
if(ldisorder_r) lkred_r=0
if(ldisorder_r2) lkred_r=0
if(ldisorder_r_bulk) lkred_r=0

lambda_par_kr=lambda_par_k!*0.1

!can set Vcf, and hybridization is different (should just be a rotation), z*2 for sbz
if(phi_type=="srb" .or. phi_type=="sbz") then
tc_=-1
tf_=-tc_
ec_=-ef_
!Vcf_=tc_	!<0
lscfloop=0
lscfloop_kr=0
lscfk=0
lscfkr=0
b_init=1
mu_init=0
lambda_init=0
endif

if(phi_type=="src") then
tc_=-1
!tf_=-tc_
ec_=0 !
!Vcf_=tc_	!<0
lscfloop=0
lscfloop_kr=0
lscfk=0
lscfkr=0
b_init=1
mu_init=0
lambda_init=0
print *, phi_type
endif

!
dist_max=sqrt(real(dist))+0.01

!
lqpi=1		!does qpi (and prints dos(r))
!ljdos=1
lphaseshifts=1	!does phaseshift
larpes=1	!prints arpes
!
lram=0		!try to save as much ram as possible for exp_factor
ldiag_g0k=1	!use eigenvectors to build G0k
lfft=1		!use fft
!
!if (phi_type=="dze") then
!lkred=.false.
!lrednscf=.false.
!else
!lkred=.true.		!reduced set of k points in scf calc
!lrednscf=.true.	!reduced set of k points in nscf calc for kr
!endif
!lkred=1
!
format_smear="(20E15.6)"
format_wr="(2i8,2e13.5)"
format_er="(i8,e15.6)"
!
!!!
!Creates directories
call system('mkdir ' //trim(label))
if (lspectrfunct .or. lspectrfunct_kr) call system('mkdir ' //trim(label)//"/sf")
if (lwritewfc)    call system('mkdir ' //trim(label)//"/wfc")
call system('cp par.in ' //trim(label))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Writes or reads parametrs

INQUIRE(file=trim(label)//'/param', EXIST=lfileexists)
OPEN(unit=32,file=trim(label)//'/param',status='unknown')
OPEN(unit=24,file=trim(label)//'/param_info',status='unknown')


if (lreadparam .and. lfileexists) then
 write(*, "(a40)") "Reading parametrs from file 'param'"
 read (32, "(5i7)") nx,ny,nz
 read (32, "(10f10.5)") ef_, ec_, tf_, tc_, Vcf_, T_min, ef0, ec0, filling
 read (32, "(10i5)") pbcx, pbcy, pbcz, lfixmu_kr, Vtki, lfixmu
 read (32, "(a3)") phi_type
 read (32, "(f10.5)") eta
! read (32, "(5i7)") nkx_r, nky_r,nkx_kk, nky_kk
! read (32, "(6i5)") lkshifted
!read (32, "(10f10.5)"), T_g
! read (32, "(10i5)"), lmagimp

 write (*, "(5i7)") nx,ny,nz
 write (*, "(10f10.5)") ef_, ec_, tf_, tc_, Vcf_, T_min, ef0, ec0, filling
 write (*, "(10i5)") pbcx, pbcy, pbcz, lfixmu_kr, Vtki, lfixmu
 write (*, "(a3)") phi_type

else
 write(*, "(a40)") "Writing parametrs to file 'param'"
 write (32, "(5i7)") nx,ny,nz
 write (32, "(10f10.5)"), ef_, ec_, tf_, tc_, Vcf_, T_min, ef0, ec0, filling
 write (32, "(10i5)"), pbcx, pbcy, pbcz, lfixmu_kr,Vtki, lfixmu
 write (32, "(a3)"), phi_type
 write (32, "(f10.5)") eta
 write (32, "(5i7)") nkx_r, nky_r, nkx_kk, nky_kk
 write (32, "(6i5)") lkshifted
 write (32, "(10f10.5)"), T_g
 write (32, "(10i5)"), lmagimp
 
 call write_info

endif

 close(unit=32)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
nkz_plot=nz
!
nkx_k=nx*nkx_kk*nkx_r
nky_k=ny*nky_kk*nky_r
!nkx_k=nx*nkx_kk
!nky_k=ny*nky_kk
nkx=nx*nkx_r
nky=ny*nky_r
nkx_kr=nkx_r/nx_kr*nx
nky_kr=nky_r/ny_kr*ny
nkz=nz
print *, "nkx_k=", nkx_k
print *, "nky_k=", nky_k
print *, "nkx=", nkx
print *, "nky=", nky
print *, "nkx_kr=", nkx_kr
print *, "nky_kr=", nky_kr
print *, "nkx_r=", nkx_r
print *, "nky_r=", nky_r
print *, "nkz=", nkz


!
if (nkx>2) then
 lkx=.true.	
 else
 lkx=.false.	
endif
if (nky>2) then
 lky=.true.	
 else
 lky=.false.	
endif
if (nkz>2) then
 lkz=.true.	
 else
 lkz=.false.	
endif
if (nkx>3) then
 lkx2=.true.	
 else
 lkx2=.false.	
endif
if (nky>3) then
 lky2=.true.	
 else
 lky2=.false.	
endif
if (nkz>3) then
 lkz2=.true.	
 else
 lkz2=.false.	
endif
if (nkx>4) then
 lkx3=.true.	
 else
 lkx3=.false.	
endif
if (nky>4) then
 lky3=.true.	
 else
 lky3=.false.	
endif
if (nkz>4) then
 lkz3=.true.	
 else
 lkz3=.false.	
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
total_size=2*nx*ny*nz*(l_orbital+l_orbital_f)
!k_total_size=4*nkx*nky*nkz
n_sites=nx*ny*nz
n_k=nkx_k*nky_k*nkz
N_el=2*n_sites*filling*l_orbital		!total electrons to be kept
N_el_k=2*n_k*filling*l_orbital		!total electrons to be kept
N_el_kr=2*nz*filling*nkx*nky*l_orbital
thr_mu=thr_mu*N_el_k
nc_0=2*filling*l_orbital-nf_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

time_init=0
time_k=0
time_buildh=0
time_diag=0
time_nf=0
time_mu=0
time_buildh_k=0
time_diag_k=0
time_wr=0
time_wrk=0
time_nfk=0
time_zhpev=0
time_kr_setvar=0
time_kr_buildh=0
time_kr_diag=0
time_kr_nf=0
time_kr_mu=0
time_kr_write=0
time_invert=0
time_matmul=0
time_time_ev=0
time_r8f2=0
time_nfimag=0
time_e=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do ix=0,nx-1
!do iy=0,nx-1
!do iz=0,nx-1
!do icf=0,1
!do is=0,1
!do o1=0,l_orb(icf)-1
!ind=index_xyzcso(ix,iy,iz,icf,is,o1,nz)
!enddo
!enddo
!enddo
!enddo
!enddo
!enddo


time_0=secnds(0.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start calculations

call open_files

if (lreadham) call readham 

call compute_nn(n_nn, n_nn2, n_nn3)		!just gives the number, to allocate variables

call allocate_var_doublet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set Kramer's doublet

 call set_cf_matrix(alpha_jzi_mat)
 call set_doublet!(alpha_jzi, alpha_msigma)
 call set_spin 

!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Sets nearest neighbours, with position wrt to zero and angles 

if (.not. lreadham) then

 call set_nn(nn_theta_phi, nn_xyz, nn_ijk)
 call set_nn2(nn2_theta_phi, nn2_xyz, nn2_ijk)
 call set_nn3(nn3_theta_phi, nn3_xyz, nn3_ijk)

 call set_phi !(phi_sigma_alpha, phi_alpha_sigma)
 call set_kin_en!(kin_energy)

 call write_nn_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)
 call write_nn2_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)
 call write_nn3_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)

endif	!lreadham

 call set_zvalues(nz_plot, z_values_plot)

time_1=secnds(0.0)

!write(*, "(a25,f10.2)") "Initialization: seconds ", time_1-time_0
!time_init=time_init+ time_1-time_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! K loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (loopk) then
time_0=secnds(0.0)

 call do_k_loop

time_1=secnds(0.0)
time_k=time_1-time_0
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kr loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (loopkr) then
time_0=secnds(0.0)

 call  do_kr_loop

time_1=secnds(0.0)
time_kr=time_1-time_0
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Real space loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (loopr) then
time_0=secnds(0.0)


!scf
if (lscfr) call do_r_loop (.true.)

!nscf
if (lnscfr) call do_r_loop (.false.)


time_1=secnds(0.0)
time_r=time_1-time_0
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call deallocate_var_doublet
call close_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(a21,f10.2)") ""

if (loopk) then
write(*,"(a21,f10.2)") "Time k:", time_k
write(*,"(a21,f10.2)") "Time set var k:", time_setvark
write(*,"(a21,f10.2)") "Time build h k:", time_buildh_k
write(*,"(a21,f10.2)") "Time diag h k:", time_diag_k
write(*,"(a21,f10.2)") "Time write to file k:", time_wrk
write(*,"(a21,f10.2)") "Time nf,mu,scf k:", time_nfk
write(*,"(a21,f10.2)") "Time nscfk:", time_nscfk
write(*,"(a21,f10.2)") ""
endif

if (loopkr) then
write(*,"(a21,f10.2)") "Time kr:", time_kr
write(*,"(a21,f10.2)") "Time kr setvar:", time_kr_setvar
write(*,"(a21,f10.2)") "Time build h kr:", time_kr_buildh
write(*,"(a21,f10.2)") "Time diag h kr:", time_kr_diag
write(*,"(a21,f10.2)") "Time mu kr:", time_kr_mu
write(*,"(a21,f10.2)") "Time nf kr:", time_kr_nf
write(*,"(a21,f10.2)") "Time write kr:", time_kr_write
write(*,"(a21,f10.2)") "Time nscfkr:", time_nscfkr
!write(*,"(a21,f10.2)") "Time rnscf:", time_rnscf
write(*,"(a21,f10.2)") ""
endif

if (loopr) then
write(*,"(a21,f10.2)") "Time r:", time_r
write(*,"(a21,f10.2)") "Time buildh:", time_buildh
write(*,"(a21,f10.2)") "Time diag:", time_diag
write(*,"(a21,f10.2)") "Time nf:", time_nf
write(*,"(a21,f10.2)") "Time mu:", time_mu
write(*,"(a21,f10.2)") "Time self-cons:", time_sc
write(*,"(a21,f10.2)") "Time write to file:", time_wr
write(*,"(a21,f10.2)") "Time spectral fcts:", time_sfr
write(*,"(a21,f10.2)") "Time Green:", time_green
write(*,"(a21,f10.2)") "Time Time_ev:", time_time_ev
write(*,"(a21,f10.2)") "Time R8f2:", time_r8f2
write(*,"(a21,f10.2)") "Time nfimag:", time_nfimag
write(*,"(a21,f10.2)") "Time compute_E:", time_E
endif
write(*,"(a21,f10.2)") ""
write(*,"(a21,f10.2)") "Time zhpev:", time_zhpev
write(*,"(a21,f10.2)") "Time invert:", time_invert
write(*,"(a21,f10.2)") "Time matmul:", time_matmul


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine do_r_loop(lscfr)
real(kind=idp)	:: time_0,time_1
logical, intent(in)	::lscfr

write (*, "(a20)") "Real space loop"

call allocate_var_r

if (n_sites<10) call write_states

T=T_min			!no T loop
one=1.0d0
zero=0.0d0


call set_siteenergies_r !sets on-site energies (for disorder)

call set_initialvalues_r	!sets initial b, lambda, mu


diff=1
i_iter=0
ok=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set k points

!scf run
if (lscfr) then			
nk_r_tot=nkx_r*nky_r
 if (lkred_r) then
  call set_kvecs_kr_pw (k_vecs_r,nk_r, nkx_r, nky_r,1, lkshifted_r)		!reduced set
 else
if(.not. ltr)  call set_kvecs2d     (k_vecs_r,nk_r, nkx_r, nky_r,.not. lkshifted_R)	!whole set, true= include gamma
if(ltr)        call set_kvecs2d_tr  (k_vecs_r,nk_r, nkx_r, nky_r,.not. lkshifted_R)	!whole set reduced by 2, true= include gamma
 endif
!nscf run, not path
else if (.not. lpathnscfr) then		
nk_r_tot=nkx_r_nscf*nky_r_nscf
 if (lkred_r) then
  call set_kvecs_kr_pw (k_vecs_r,nk_r, nkx_r_nscf, nky_r_nscf,3,lkshifted_r_nscf)		!reduced set
 else
if(.not. ltr)  call set_kvecs2d     (k_vecs_r,nk_r, nkx_r_nscf, nky_r_nscf,.not. lkshifted_R_nscf)	!whole set, true= include gamma
if(ltr)        call set_kvecs2d_tr  (k_vecs_r,nk_r, nkx_r_nscf, nky_r_nscf,.not. lkshifted_R_nscf)	!whole set reduced by 2, true= include gamma
 endif
 if (mod(nkx_r_nscf,2)==1) call set_kvecs_r_arpes(k_vecs_r_arpes,nk_r_arpes,list_k_vecs_r_arpes)	!nscf for arpesy if odd number of k points
!nscf run, path
 else					
 call set_kr_path_r(k_vecs_r,nkpath_r,nk_r)	!nscf for k path
endif

write (*, "(a10)"), "~~~~~~~"

write (*, "(a10,i5)"), "nk_r=", nk_r, "n_k_kr=",n_k_kr

!do ik=0, nk_r-1
!write (*, "(4f15.8)"), k_vecs_r(ik,1), k_vecs_r(ik,2), k_vecs_r(ik,3), k_vecs_r(ik,4)
!enddo

write (*, "(a10)"), "~~~~~~~"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Builds V (needed in SCF step) and allocates wfc_tot

lwfctot=.false.
if (lscfr) lwfctot=.true.		!SCF
if (.not. lscfr .and. .not. lpathnscfr .and. lrhoktot) lwfctot=.true.	!NSCF, keeping al k in memory
if (.not. lscfr .and.       lpathnscfr .and. lsfrkred) lwfctot=.true.	!PATH, reduced procedure
if (.not. lscfr .and.       lpathnscfr .and. .not. lfullbz	) lwfctot=.true.	!PATH, reduced brilloin zone



if (lwfctot) call build_V_r(h_psi_tot, index_psi_tot)			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Main self consistent loop


do while (diff>threshold .and. (lscfloop .or. .not. lscfr))

  time_0=secnds(0.0)

  i_iter=i_iter+1

if (lscfr)  write(*,"(i5,a10)") i_iter, "SCF loop"
if (.not. lscfr)  write(*,"(a10)") "NSCF"

!!!

if ( lwfctot )  e_r_tot=0
if ( lwfctot )  w_r_tot=0

! k loop

  do ik=0, nk_r-1

  write (*,"(a,i5,4f10.5)") "k point", ik ,k_vecs_r(ik,1),k_vecs_r(ik,2),k_vecs_r(ik,3),k_vecs_r(ik,4)

!!!!!!!!!!!!!!!!Allocate couplings

!call build_couplings(b,ec,ef,tc,tf,Vcf)		!uso b_old

!!!!!!!!!!!Allocate Hamiltonian

!print *, "build_ham"
  if (.not. lreadwfc .or. i_iter>1) call build_ham_r(ik,Ham,b,lambda,mu,ec_site,ef_site)			!uso b_old, lambda_old, mu_old
 
  if (n_sites < 8.and. i_iter==1) call write_Ham(ham, "Ham")

!!!!!!!!!!Rotates into k space, for check purposes
 
 ! if (n_sites < 5) then 
  
!  call set_kvecs (k_vecs)  
  
 !  call set_U(U,Um1)		!matrix for rotation to k space
!
  ! if (i_iter==1) call write_Ham(U,"U")
 
   !call BASIS_ROTATIONcmplx (U, Ham, Um1, Ham_rot)
  
 !  if (i_iter==1) call write_Ham(Ham_rot, "Ham_rot")

!  endif


time_1=secnds(0.0)
time_buildh=time_buildh+time_1-time_0
 
!!!!!!!!!!!Diagonalize Hamiltonian

time_0=secnds(0.0)

!print *, "diag_ham"
 if (.not. lreadwfc .or. i_iter>1) call diagonalize_ham(Ham,e,w)
 if (lreadwfc .and. i_iter==1) call readwr(ik,e,w)

time_1=secnds(0.0)
time_diag=time_diag+time_1-time_0

  if (n_sites < 8 .and. i_iter==1) call write_energies (e,w)

!stop

!print *, lscfr,.not. lpathnscfr,lsfrkred,.not. lfullbz

if ( lwfctot )  w_r_tot(:,:,ik)=w(:,:)
if ( lwfctot )  e_r_tot(:,ik)=e(:)



!PATH, no reduced procedure (one K at a time)
!if (.not. lscfr .and. lpathnscfr .and. .not. lsfrkred )  e=e-mu
if (.not. lscfr .and. lpathnscfr .and. .not. lsfrkred )  call spectralfunction_r_xyzk_zk(ik, lpathnscfr)
!if (.not. lscfr .and. lpathnscfr .and. .not. lsfrkred )  call spectralfunction_r_k(ik)

!NSCF, one K at a time
!if (.not. lscfr .and. .not. lpathnscfr .and. .not. lrhoktot )  e=e-mu
if (.not. lscfr .and. .not. lpathnscfr .and. .not. lrhoktot )  call spectralfunction_r_xyzk_zk(ik, lpathnscfr)

!print *, "a"


  enddo	!kpoints

!!!!!!!!!!!!Compute chemical potential mu

time_0=secnds(0.0)

  if (lfixmu .or. .not. lscfr) then
   mu_new=mu
  else  
   call compute_mu_r(e_r_tot)									!da fare con k points
  endif

time_1=secnds(0.0)
time_mu=time_mu+time_1-time_0

!!!!!!!!Compute occupations, new b/lambda, and total free energy

time_0=secnds(0.0)
 
 if (lscfr .or. .not. lpathnscfr .and. lwfctot)  &
 call compute_nf_r (nf,nc,nfc,nff,b_new,lambda_new,b,lambda,0.0d0,e_r_tot,w_r_tot,diff, k_vecs_r,nk_r,nk_r_tot, total_size,n_sites,nx,ny,nz,h_psi_tot,index_psi_tot,lvac_r, lambda_par)
! if (lscfr .or. .not. lpathnscfr .and. lwfctot)  call compute_nf_r (nf,nc,b_new,e_r_tot,w_r_tot,ef, lambda_new, nfc,nff,diff,0.0d0)		!ottengo b_new

time_1=secnds(0.0)
time_nf=time_nf+time_1-time_0


!!!!!!!! Prints results to file
!scf b , lambda, mu

if (lscfr) then
 OPEN(unit=37,file=trim(label)//'/b_r2',status='unknown')

 write(37,"(3i4)") nx,ny,nz

 do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z

 ind=ix + iy*nx + iz*nx*ny		
 write(37,"(3i4,100f15.10)") ix,iy,iz, b_new(ind), lambda_new(ind)

   enddo
  enddo
 enddo
 write(37,"(f20.10)")  mu_new

 close (unit=37)
endif

!!!!!!!!!!!Check selfconsistency

time_0=secnds(0.0)





  !!!!!!!!!!!!!Set new b and lambda

  !call random_seed
  !call random_number(delta_alpha)
  !alpha=(alpha_min+(alpha_max-alpha_min)*delta_alpha)!/(real(i_iter)**(0.6))
  !print *, "alpha(mu and b)=", alpha

  if (.not. lscfr) diff=0	!I just do the first iteration
  
  if (diff<threshold) then
   print *, "Done!"
   ok=.true.
  else

   alpha=alpha_min
   
   b=alpha*b_new+(1-alpha)*b
   lambda=alpha*lambda_new+(1-alpha)*lambda
   mu=alpha*mu_new+(1-alpha)*mu

  endif


time_1=secnds(0.0)
time_sc=time_sc+time_1-time_0

  print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ending main loop


!!!
!writes to file
!eigenvectors

!if (lscfr.and. (.not. lreadwfc .or. i_iter >1)) call write_wrtot (e_r_tot,w_r_tot)
time_0=secnds(0.0)
if ( lwritewfc .and. lwfctot .and. (.not. lreadwfc .or. i_iter >1)) call write_wrtot (e_r_tot,w_r_tot)
time_1=secnds(0.0)
time_wr=time_1-time_0


!!!

!call write_info

!!!
!print *, lscfr

!stop

if (ok .and. lscfr) then
! write(*,"(3a4,100a15)") "#x", "y", "z", "b","nc", "nf","ef-lambda-mu", "lambda", "nfc", "nff"
 write(11,"(3a4,100a15)") "#x", "y", "z", "b","nc", "nf","ef-lambda-mu", "lambda", "nfc", "nff"

 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
   do iz=0,nz-1	!z

 ind=ix + iy*nx + iz*nx*ny		
! write(*,"(3i4,100f15.10)") ix,iy,iz, b_new(ind),  nc(ind), nf(ind),ef(ind),lambda_new(ind),nfc(ind), nff(ind)
 write(11,"(3i4,100f15.10)") ix,iy,iz, b_new(ind),  nc(ind), nf(ind),ef(ind,0),lambda_new(ind),nfc(ind), nff(ind)

   enddo
  enddo
!  write(11,"(a1)"), " " 
 enddo

 write(11,"(a10,i10)")  "#n_iter=", i_iter



 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
   do iz=0,0	!z

 ind=ix + iy*nx + iz*nx*ny		
! write(*,"(3i4,100f15.10)") ix,iy,iz, b_new(ind),  nc(ind), nf(ind),ef(ind),lambda_new(ind),nfc(ind), nff(ind)
 write(11,"(3i4,100f15.10)") ix,iy,iz, b_new(ind),  nc(ind), nf(ind),ef(ind,0),lambda_new(ind),nfc(ind), nff(ind)

   enddo
  enddo
!  write(11,"(a1)"), " " 
 enddo



 OPEN(unit=37,file=trim(label)//'/b_r_ok',status='unknown')

 write(37,"(3i4)") nx,ny,nz

 do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z

 ind=ix + iy*nx + iz*nx*ny		
 write(37,"(3i4,100f15.10)") ix,iy,iz, b_new(ind), lambda_new(ind)

   enddo
  enddo
 enddo
 write(37,"(f20.10)")  mu_new

 close (unit=37)


 if (loopk) write(*,"(a3,6f14.8,i5)")  "k:",T,b_new_k, lambda_new_k, mu_new_k, E_tot_k

endif


!!!!!!!!!!!!!
!POSTPROCESSING, if I have wfc_tot
!!!

if (lwfctot) then

 time_0=secnds(0.0)

 call write_bands_r (e_r_tot, lscfr)		!BANDS

! call part_ratio(e_r_tot, w_r_tot, lscfr)	!PARTICIPATION RATIO

 if (lwfc) call write_wfc(e_r_tot,w_r_tot)	!WFC

 if (lspectrfunct) then
 
  !if (lscfr .or. .not. lpathnscfr) 
  if(.not. lscfr .and. .not. lpathnscfr) call spectralfunction_r(e_r_tot,w_r_tot)	!NSCF local DOS nx,ny,nz (summed over k)


if (lpathnscfr)  n_ksf=n_k_total
if (lscfr .or. .not. lpathnscfr)  n_ksf=nk_r

  do ik=0, n_ksf-1
   if(.not. lscfr) call spectralfunction_r_xyzk_zk(ik, lpathnscfr)	!NSCF, PATH, k resolved
  enddo
 
 endif


!  if (lpathnscfr  .and. lsfrkred) then
!   do ik=0, n_k_total-1
!    call spectralfunction_r_k(ik)		!PATH DOS nz,k
!  endif


 time_1=secnds(0.0)
 time_sfr=time_1-time_0

endif
!!!
!END OF POSTPROCESSING
!!!!!!!!!

 call do_dissipation

time_time_ev_0=secnds(0.0)
 call time_evolution
time_time_ev=secnds(0.0)-time_time_ev_0

!!!!!!!!!
!Green's function 

print *, "ok"
 time_green_0=secnds(0.0)

if (lscfr .and. (lgreen .or. lborn .or. lg0)) call do_green


 time_green_1=secnds(0.0)
 time_green=time_green_1-time_green_0

!!!!!

call deallocate_var_r

end subroutine do_r_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine do_k_loop
real(kind=idp)	:: time_0,time_1


write(*,"(a25)") ""
write(*,"(a25)") "K loop"
write(*,"(a25)") ""

time_0=secnds(0.0)

if (lkred) then
 call set_kvecs_pw (k_vecs3d, nk3d_irr,nkx_k,nky_k,nkz,lysameasx,lzsameasx,nky)	!reduced set
!call set_kvecs_pw (k_vecs3d, nk3d_irr,nkx_k,nkx_k,nkx_k,.true.,.true.,nkx_k)	!reduced set
! call set_kvecs_cubic(k_vecs3d,nk3d_irr,nkx_k)
nk3d=nk3d_irr
else
 call set_kvecs (k_vecs3d)				!whole set
nk3d=n_k
endif

!!!!
if (nk3d<100) then
do ik=0, nk3d-1
write(*, "(i5,1x,3f20.10,f7.2)"), ik, k_vecs3d(ik,1),k_vecs3d(ik,2),k_vecs3d(ik,3),k_vecs3d(ik,4)
enddo
endif
!!!!

 call allocate_var_k


if (.not. lreadham) call set_phi_k(phi_k_sigma_alpha_tot,phi_k_alpha_sigma_tot,kin_energy_k_tot, delta_k_tot,nk3d,k_vecs3d,3)

!if (.not. lreadham) call factor_k(factor_k_tot,nk3d,k_vecs3d,3)


if (nk3d<100) then
 call write_k_points(k_vecs3d)
endif


 call build_v_k(v_k_tot)

write(*,"(a35,f15.5, 2i10)") "N_el_k, n_k, n_k_used=", N_el_k, n_k, nk3d
write(*,"(a35,2f15.5)") "nf_0,nc_0=", nf_0,nc_0

write(*,"(3(a10,f15.8))") "T_min=", t_min,"T_max=", t_max

 call set_initialvalues_k


time_1=secnds(0.0)
time_setvark=time_1-time_0

!t=t_min

!print *, t_min, t_max, delta_t

!T loop
do T=T_min, T_max, delta_T
 write(*,"(3(a10,f15.8))") "T=", T

!only sets mu
!(1)

if (lmuscf .and. .not. lscfk .and. .not. lmufscf) then
print *, "setting mu nscf only"
!   b_k=1
!   lambda_k=0

  e_k_tot=0		!vector of all eigenvalues
  w_k_tot=0
 
  do ik=0, nk3d-1
 
   call build_ham_k(ik,Ham_k,b_k,lambda_k, 0.0d0,.true.,nk3d,phi_k_sigma_alpha_tot,phi_k_alpha_sigma_tot, kin_energy_k_tot,k_vecs3d)		!I put mu=0

   call diagonalize_ham_k(Ham_k,e_k,w_k)

   if (n_k< 10 .and. i_iter==1)   call write_energies_k (e_k,w_k)

   e_k_tot(:,ik)=e_k(:)
   w_k_tot(:,:,ik)= w_k(:,:)

  enddo		!k points

   call compute_mu_k_tot(e_k_tot,mu_k)


!
!(2)
elseif(lmuscf .and. lmufscf .and. .not. lscfk) then
print *, "setting mu and lambda nscf"

!scf-loop
 diff=1
 i_iter=0
 ok=.false.

 do while (diff>100*threshold_k)
  
  i_iter=i_iter+1
  if (i_iter >5e5) stop

  e_k_tot=0		!vector of all eigenvalues
  w_k_tot=0
 
! k-vecs loop
  do ik=0, nk3d-1
 
   time_0=secnds(0.0)
 
   call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk3d,phi_k_sigma_alpha_tot,phi_k_alpha_sigma_tot, kin_energy_k_tot,k_vecs3d)		!uso b_old, lambda_old, mu_old

  time_1=secnds(0.0)
  time_buildh_k=time_buildh_k+time_1-time_0

  time_0=secnds(0.0)

   call diagonalize_ham_k(Ham_k,e_k,w_k)

   e_k_tot(:,ik)=e_k(:)
   w_k_tot(:,:,ik)= w_k(:,:)

  time_1=secnds(0.0)
  time_diag_k=time_diag_k+time_1-time_0


  enddo		!k points

!!!!!!!!!!!!Compute chemical potential mu

  time_0=secnds(0.0)

   call compute_mu_k(e_k_tot)											!ottengo mu_new
   call compute_nfnc_k (nf_k,nc_k,b_new_k,e_k_tot,w_k_tot,ef_k, lambda_new_k, nfc_k, nff_k,diff, 0.0d0)		!ottengo mu_new and lambda_new


if (abs(alpha_max_k-alpha_min_k)<1e-5) then
alpha_k=alpha_min_k
else
 call random_seed
 call random_number(delta_alpha)
alpha_k=(alpha_min_k+(alpha_max_k-alpha_min_k)*delta_alpha)!/(real(i_iter)**(0.6))
endif

if ((T<T_min+0.01) .and. (i_iter < 100 .or. diff>diff_thr)) then

alpha_k_u= alpha_k*0.005
alpha_mu_k_u= alpha_mu_k*0.01
else
alpha_k_u=alpha_k
alpha_mu_k_u=alpha_mu_k
endif

!print *, "alpha(mu and b)=", alpha_k, alpha_k_u


  b_k=alpha_k_u*b_new_k+(1-alpha_k_u)*b_k
  lambda_k=alpha_k_u*lambda_new_k+(1-alpha_k_u)*lambda_k
  mu_k=alpha_mu_k_u*mu_new_k+(1-alpha_mu_k_u)*mu_k


 enddo	!SCF loop

!just to get the nf file
e_k_tot=e_k_tot+mu_k
   call compute_mu_k_tot(e_k_tot,mu_k)
e_k_tot=e_k_tot-mu_k


 time_0=secnds(0.0)


!!!!!!!!!!!!!!!!!!!(3)
elseif(lscfk) then

!scf-loop
 diff=1
 i_iter=0
 ok=.false.

 do while (diff>threshold_k)
  
  i_iter=i_iter+1
  if (i_iter >5e5) stop

  e_k_tot=0		!vector of all eigenvalues
  w_k_tot=0
 
! k-vecs loop
  do ik=0, nk3d-1
 
   time_0=secnds(0.0)
 
   call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk3d,phi_k_sigma_alpha_tot,phi_k_alpha_sigma_tot, kin_energy_k_tot,k_vecs3d)		!uso b_old, lambda_old, mu_old

  ! if (i_iter==1)  write(*,"(i5,a15,3(a10,f10.5))") ik, "k point", "F_k", factor_k_tot(ik),  "weight", k_vecs3d(ik,4)

 !  if (nk3d_irr< 10 .and. i_iter==1) call write_Ham_k(ham_k)

  time_1=secnds(0.0)
  time_buildh_k=time_buildh_k+time_1-time_0

  time_0=secnds(0.0)

   call diagonalize_ham_k(Ham_k,e_k,w_k)

   if (n_k< 10 .and. i_iter==1)   call write_energies_k (e_k,w_k)

   e_k_tot(:,ik)=e_k(:)
   w_k_tot(:,:,ik)= w_k(:,:)

  time_1=secnds(0.0)
  time_diag_k=time_diag_k+time_1-time_0

  enddo		!k points

!!!!!!!!!!!!Compute chemical potential mu

  time_0=secnds(0.0)

  if (lfixmu_k) then
  print *, "fixed mu"
   mu_new_k=mu_k
  else  
  print *, "calc mu"
   call compute_mu_k(e_k_tot)
  endif

!print *, mu_k, mu_new_k

!!!!!!!!Compute occupations, new b/lambda, and total free energy and Check selfconsistency

  call compute_nf_k (nf_k,nc_k,b_new_K,e_k_tot,w_k_tot,ef_k, lambda_new_k, nfc_k, nff_k,diff, 0.0d0)		!ottengo b_new
!  call compute_nf_k (nf_k,nc_k,b_new_K,e_k_tot,w_k_tot,ef_k, lambda_new_k, nfc_k, nff_k,diff, mu_new_k)	!ottengo b_new

!!!!!!!!!!!Check selfconsistency

!  diff=(1-b_new_k**2-nf_k)**2+((Nfctot-N_el_k)/n_k)**2
!if ( .not. lfixmu_k .and. lmucf)  	diff=(1-b_new_k**2-nf_k)**2+((Nfctot-N_el_k))**2+(b_new_k-nfc_k/2/(lambda_new_k-nff_k))**2
!if ( .not. lfixmu_k .and. .not. lmucf) 	diff=(1-b_new_k**2-nf_k)**2+((Nctot-N_el_k/2))**2+(b_new_k-nfc_k/2/(lambda_new_k-nff_k))**2
!if (lfixmu_k)  		diff=(1-b_new_k**2-nf_k)**2+(b_new_k-nfc_k/2/(lambda_new_k-nff_k))**2

!  diff=sqrt(diff)
!!!!!!!!!!!!!Set new b and lambda


  if (diff<threshold_k) then
  ok=.true.
  else
!if (i_iter < 1000 .and. T==T_min) then

if (abs(alpha_max_k-alpha_min_k)<1e-5) then
alpha_k=alpha_min_k
else
call random_seed
call random_number(delta_alpha)
alpha_k=(alpha_min_k+(alpha_max_k-alpha_min_k)*delta_alpha)!/(real(i_iter)**(0.6))
endif

if ((T<T_min+0.01) .and. (i_iter < 10 .or. diff>diff_thr)) then

alpha_k_u= alpha_k*0.005
alpha_mu_k_u= alpha_mu_k*0.01
else
alpha_k_u=alpha_k
alpha_mu_k_u=alpha_mu_k
endif

!print *, "alpha(mu and b)=", alpha_k, alpha_k_u


  b_k=alpha_k_u*b_new_k+(1-alpha_k_u)*b_k
  lambda_k=alpha_k_u*lambda_new_k+(1-alpha_k_u)*lambda_k
  mu_k=alpha_mu_k_u*mu_new_k+(1-alpha_mu_k_u)*mu_k
  endif

  time_1=secnds(0.0)
  time_nfk=time_nfk+time_1-time_0

 enddo	!SCF loop



 time_0=secnds(0.0)

endif	!lscfk

! print *, phi_type

 if ((lscfk .and. abs(T-T_min)<1e-5) .or. phi_type=="sra".or. phi_type=="srb".or. phi_type=="sbz".or. phi_type=="sxy" .or. lmuscf) then
! print *, phi_type
  OPEN(unit=30,file=trim(label)//'/b_k2',status='unknown')
  write(30,"(3f20.10)")  b_k, lambda_k,  mu_k	!only if I did a scf run
  close (unit=30)
 endif


print *, "Set kvecs nscf"
 if (lnscfk) 	 	     call set_k_vecs_nscf!(k_vecs_nscf,factor_k_nscf,phi_k_sigma_alpha_nscf,phi_k_alpha_sigma_nscf,delta_k_nscf)
 if (lnscfk) 		     call set_k_vecs_path3d!
 if (lnscfk) 		     call set_k_vecs_hsp!
! if (lnscfk) 		     call set_k_vecs_chern



!stop

 if (lnscfk) call do_k_nscf(b_k,lambda_k, mu_k)

 time_1=secnds(0.0)
 time_nscfk=time_1-time_0

 time_0=secnds(0.0)

 if (ok) then

 write(*,"(a20)")  " "
if (t==t_min) write(13,"(14a11,a10)")  "T","b_k","bnfc", "nf_k" , "ef_k", "b**2+nf","nfc_k", "lambda_k",   "nc_k","nf_k+nc_k", "nff_k","mu_k","nf7_k","i_iter"
 write(13,"(13f11.6,i10)")  T, b_k,(nfc_k/(lambda_k-nff_k)/2), nf_k,  ef_k, b_new_k**2+nf_k, nfc_k, lambda_k,  nc_k, nf_k+nc_k, nff_k, mu_k, nf7_k, i_iter

!print *, nctot, n_el_k/2, nfctot, n_el_k
 
 e_k_tot=e_k_tot-mu_new_k!+mu_k


if (abs(T-T_min)<1e-6) call write_bands_k(e_k_tot)

if (abs(T-T_min)<1e-6 .and. lspectrfunct_kr) call spectralfunction_k(e_k_tot,w_k_tot)


 endif

time_1=secnds(0.0)
time_wrk=time_1-time_0

enddo	!T loop

call deallocate_var_k

write(*,"(a25)") ""
write(*,"(a25)") "end K loop"
write(*,"(a25)") ""



end subroutine do_k_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_info
integer	::i

!OPEN(unit=24,file=trim(label)//'/param_info',status='unknown')

do i=6,24,18

write (i, "(a60,f15.5,4i5,f10.5)"), "#N_el, size (nx,ny,nz,ntot), filling (1 is half filling)=", N_el, nx,ny,nz, n_sites,real(N_el)/real(n_sites)/2.
write (i, "(a4,f10.5)"), "#Ef=", ef_
write (i, "(a4,f10.5)"), "#Ec=", ec_
write (i, "(a4,f10.5)"), "#tf=", tf_
write (i, "(a4,f10.5)"), "#tc=", tc_
write (i, "(a4,f10.5)"), "#Vcf=", Vcf_
write (i, "(a4,f10.5)"), "#T=", T_min
write (i, "(a6,i5)"), "#Vtki=", Vtki
write (i, "(a4,f10.5)"), "#mu=", mu
write (i,"(a10,100f15.10)") "#<nf>,<nc>= ", Nftot/n_sites, Nctot/n_sites
write (i,"(a20,100f15.10)") "#ntot(site, total)= ", Nfctot/n_sites, Nfctot
write (i, "(a4,f15.10)"), "#T_K~", T_K
write (i, "(a7,i6)"), "#N_iter=", i_iter
write (i,"(a10,100f20.10)") "#!energy=", E_tot, E_tot_l

enddo

! close(unit=24)

end subroutine write_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_V_r(h_psi_tot, index_psi_tot)
integer			::index_psi_tot(0:total_size-1,n_max_link,3)
complex(kind=idpc),allocatable	::h_psi_tot(:,:,:)

print *, "build_v_r"
print *, "total_size=", total_size

allocate(V(0:total_size-1,0:total_size-1))
allocate(h_psi_tot(0:total_size-1,n_max_link,0:nk_r-1))
allocate(w_r_tot(0:total_size-1,0:total_size-1,0:nk_r-1))
allocate(e_r_tot(0:total_size-1,0:nk_r-1))


do ik=0, nk_r-1

 V=0
 
 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
   do iz=0,nz-1	!z
    do icf=0,1	!c/f
     do is=0,1	!up/down 
      do o1=0,l_orb(icf)-1
      
      !indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz + o1*4*nx*ny*nz
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx, ny, nz)

      call ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,one,zero,0.0d0,.false.,nx,ny,nz,ec_site,ef_site, k_vecs_r,nk_r)

      do i=1,n_max_link
       ind2p=index_psi(i,7)
       V(indp,ind2p)=V(indp,ind2p)+h_psi(i)
       h_psi_tot(indp,i,ik)=h_psi(i)
       index_psi_tot(indp,i,1)=index_psi(i,6)	!ind2	xyz
       index_psi_tot(indp,i,2)=index_psi(i,7)	!ind2p	xyzcso
       index_psi_tot(indp,i,3)=index_psi(i,4)	!icf	a
      enddo
      
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

 do i=0, total_size-1
 do j=i, total_size-1
 if (abs(V(j,i)-conjg(V(i,j)))>1e-9) then		! v is only for check reasons
 print *, "error in V", ik,i, j, V(i,j), V(j,i)
 !stop
 endif
 enddo
 enddo

!print *, V
enddo	!k


deallocate(v)

print *, "ok build_v_r"

end subroutine build_V_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_ham_r(ik,Ham,b, lambda,mu,ec_site,ef_site)
complex(KIND=idpc), intent(out) 	:: Ham(0:total_size-1,0:total_size-1)
real(KIND=idp),intent(in) 	:: b(0:n_sites-1), lambda(0:n_sites-1),mu,ef_site(0:n_sites-1,0:l_orbital_f-1),ec_site(0:n_sites-1,0:l_orbital-1)
integer, intent(in)::ik



ham=0

do ix=0,nx-1	!x
 do iy=0,ny-1	!y
  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1
     
    ! indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz + o1*4*nx*ny*nz
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx, ny, nz)
     
      call ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,b,lambda,mu,.true.,nx,ny,nz,ec_site,ef_site, k_vecs_r,nk_r)
!     call ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,one,zero,mu,.false.)


!write (*, "(7i4)") ik,ix,iy,iz,icf,is,indp

     do i=1,n_max_link
      ind2p=index_psi(i,7)
      Ham(indp,ind2p)=Ham(indp,ind2p)+h_psi(i)
!write (*, "(2i4,2f10.5)") i,ind2p,real(h_psi(i)), aimag(h_psi(i))
     enddo

     enddo
    enddo
   enddo
  enddo
 enddo
enddo




do i=0, total_size-1
do j=i, total_size-1
if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9) then
print *, "error in Hamiltonian", i, j, Ham(i,j), Ham(j,i)
!stop
endif
enddo
enddo

end subroutine build_ham_r





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_ham(Ham, string)
 character(len=*), intent(in)	::string
#ifdef cmplx
complex(KIND=idpc) 	:: Ham(0:total_size-1,0:total_size-1)
#else
real(KIND=idp) 	:: Ham(0:total_size-1,0:total_size-1)
#endif

write(*,"(a10)"),string

!do, i=0,total_size-1
!    write(*,"(100(a1,f5.2,a1,f5.2,a1))") ( "(",real(Ham(i,j)),",",aimag(Ham(i,j)),")", j=0,total_size-1 )
!enddo
 write(*,"(100i5)") ( j, j=-1,total_size-1 )
do, i=0,total_size-1
!    write(*,"(100f7.4)") ( abs(Ham(i,j)), j=0,total_size-1 )
    write(*,"(i5,100f5.2)") i,( sqrt(real(Ham(i,j))**2+aimag(Ham(i,j))**2), j=0,total_size-1 )
enddo

    write(*,"(a1)") ""


end subroutine write_ham


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagonalize_ham(Ham, E, W)
#ifdef cmplx
complex(KIND=idpc) 	:: Ham(0:total_size-1,0:total_size-1),W(0:total_size-1,0:total_size-1)
#else
real(KIND=idpc) 	:: Ham(0:total_size-1,0:total_size-1),W(0:total_size-1,0:total_size-1)
#endif
real(KIND=idp)		:: E(0:total_size-1)

CALL DIAGONALIZE( 'V', total_size, total_size, Ham, E, W )


end subroutine diagonalize_ham

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DIAGONALIZE( JOBZ, N, M, H, E, W )

IMPLICIT   NONE

INTEGER,          PARAMETER     :: IL = 1       
REAL(KIND=idp),   PARAMETER     :: ABSTOL = 0.0D0
 CHARACTER(LEN=1), PARAMETER     :: UPLO  = 'L'  
 CHARACTER(LEN=1), PARAMETER     :: RANGE = 'I'   
                                                
                                               
 CHARACTER, INTENT(IN)           :: JOBZ*1
INTEGER,   INTENT(IN)           :: N, M
#ifdef cmplx
COMPLEX(KIND=idpc), INTENT(IN)      :: H(N,N)
COMPLEX(KIND=idpc), INTENT(OUT)     ::  W(N,M)
#else
real(KIND=idp), INTENT(IN)      :: H(N,N)
real(KIND=idp), INTENT(OUT)     ::  W(N,M)
#endif

REAL(KIND=idp), INTENT(OUT)     :: E(N)

INTEGER                         :: M_OUT, INFO, IU, i, j
REAL(KIND=idp)                  :: VL, VU       ! Not referenced for RANGE=I

INTEGER,        DIMENSION(:),   ALLOCATABLE     :: IWORK, IFAIL
#ifdef cmplx
COMPLEX(KIND=idpc), DIMENSION(:),   ALLOCATABLE :: WORK,AP
#else
REAL(KIND=idp), DIMENSION(:),   ALLOCATABLE :: WORK,AP
#endif	
REAL(KIND=idp), DIMENSION(:),   ALLOCATABLE     :: RWORK



ALLOCATE( IWORK(5*N)     )
ALLOCATE( RWORK(7*N)     )
ALLOCATE( IFAIL(N  )     )
ALLOCATE(  WORK(8*N)     )	
ALLOCATE(  AP(N*(N+1)/2) )




IU = M          ! Maximum eigenvalue calculated

DO i=1,N
    DO j=1,i
        AP(i + (j-1)*(2*N-j)/2) = H(i,j)
    ENDDO
ENDDO

#ifdef cmplx
 !CALL ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, &
 !       ABSTOL, M_OUT, E, W, N, WORK,RWORK,IWORK, IFAIL, INFO )
 time0_zhpev=secnds(0.0)
 CALL ZHPEV( JOBZ,UPLO, N, AP, E, W, N, WORK,RWORK, INFO )
 time1_zhpev=secnds(0.0)
 time_zhpev=time_zhpev+time1_zhpev-time0_zhpev
 
#else
!SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
call dspev(jobz, uplo,n,ap,e,w,n,work,info)

!CALL DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, &
!        ABSTOL, M_OUT, E, W, N, WORK, IWORK, IFAIL, INFO )
#endif	

!IF (INFO.ne.0) THEN
!    WRITE(6,*) ' WARNING: ZHPEVX INFO = ', INFO
!    DO i=1,M
!        WRITE(6,*) '          i, IFAIL(i) = ', i, IFAIL(i)
!    ENDDO
!ENDIF


deallocate (IWORK)
deallocate (rWORK)
deallocate (Ifail)
deallocate (WORK)
deallocate (ap)

END SUBROUTINE DIAGONALIZE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE DIAGONALIZE_s(N, S, E )

IMPLICIT NONE

INTEGER,          PARAMETER     :: IL = 1       
REAL(KIND=idp),   PARAMETER     :: ABSTOL = 0.0D0
                                                
                                               
! CHARACTER, INTENT(IN)           :: JOBZ*1
INTEGER,   INTENT(IN)           :: N
COMPLEX(KIND=idpc), INTENT(IN)      :: S(N,N)
!COMPLEX(KIND=idpc), INTENT(OUT)     ::  W(N,N)
!COMPLEX(KIND=idpc), INTENT(IN)      :: vl(1,1)

complex(KIND=idpc), INTENT(OUT)     :: E(N)

INTEGER                         :: M_OUT, INFO, IU, i, j
!REAL(KIND=idp)                  :: VL, VU       ! Not referenced for RANGE=I

complex(KIND=idpc),allocatable:: copy(:,:)

INTEGER,        DIMENSION(:),   ALLOCATABLE     :: IWORK, IFAIL
COMPLEX(KIND=idpc), DIMENSION(:),   ALLOCATABLE :: WORK,AP, vl(:,:)
REAL(KIND=idp), DIMENSION(:),   ALLOCATABLE     :: RWORK
integer:: lwork

lwork=4*N

ALLOCATE( RWORK(2*N)     )
ALLOCATE(  WORK(lwork)     )	
allocate(copy(n,n))

copy=S


  allocate( vl(N,N) )


! time0_zhpev=secnds(0.0)
  call zgeev('N','N', N, copy, N, E, vl, N, vl, &
              N, work, lwork, rwork, info)
              
!pwcond              
!  lwork = 4*nchanl
!  allocate( work(lwork) )
!  allocate( rwork(lwork) )
!  call zgeev('N','V', 2*nchanl, smat, 2*nchanl, exphi, vl, 2*nchanl, vl, &
!              2*nchanl, work, lwork, rwork, info)              
              
! time1_zhpev=secnds(0.0)
! time_zhpev=time_zhpev+time1_zhpev-time0_zhpev



deallocate (rWORK)
deallocate (WORK)
deallocate (vl)
deallocate(copy)

END SUBROUTINE DIAGONALIZE_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_energies(e,w)
real(KIND=idp) 	:: E(0:total_size-1)
#ifdef cmplx
complex(KIND=idpc) 	:: w(0:total_size-1,0:total_size-1)
#else
real(KIND=idp) 	:: w(0:total_size-1,0:total_size-1)
#endif

do, i=0,total_size-1
    write(*,"(100f15.10)") ( e(i) )
enddo

!do, i=0,total_size-1
!   write(*,"(100f10.5)") ( real(w(i,j)), j=0,total_size-1 )
!enddo



END SUBROUTINE write_energies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine compute_mu_r (e_r_tot)
real(KIND=idp),intent(in) 		:: E_r_tot(0:total_size-1,0:nk_r-1)
real(KIND=idp)::alpha_mu_r

mu_old=mu

!if (i_iter==1) then 
!mu=(e_r_tot(N_el-1,0)+e_r_tot(N_el,0))/2
!endif

i_iter_mu=0
diff_mu=100
 
 
write(*,"(a40,100f15.10)") "******* Mu loop, starting mu=", mu
write(*,"(a5,5(a15),a20,a15)") "i_mu","mu",  "Nfctot",  "N_el", "diff_mu", "mu_new", "dndmu"

!!!!!!!mu loop

do while (diff_mu>thr_mu)
i_iter_mu=i_iter_mu+1

  nfctot=0
  do ik=0, nk_r-1
   do i=0, total_size-1	!eigenvalues
    nfctot=nfctot+1/(exp((e_r_tot(i,ik)-mu)/T)+1)
   enddo
  enddo
!!!!!!Difference between wanted and obtained N_el

diff_mu=abs(real(N_el)-Nfctot)

!!!!!!!!!!!!! Change mu

dndmu=0
do ik=0, nk_r-1
 do i= 0,total_size-1
  dndmu=dndmu+1/T/(exp((e_r_tot(i,ik)-mu)/T)+1)/(exp((-e_r_tot(i,ik)+mu)/T)+1)		!d(n)/d(mu)=sum_i f'(e_i-mu)
 enddo
enddo

dndmu=max(dndmu,0.01)		!to avoid dndmu --> 0 when T-->0


mu_new=mu+(N_el-Nfctot)/dndmu

write(*,"(i5,5f15.10,f20.10, f15.10)") i_iter_mu, mu,  Nfctot,  real(N_el), diff_mu, mu_new, dndmu


	if (diff_mu>0.5) then
!	thresh=0.1
	alpha_mu_r=alpha_mu*0.01
	else
!	thresh=thr_mu
	alpha_mu_r=alpha_mu
	endif


mu_new=alpha_mu_r*mu_new+(1-alpha_mu_r)*mu
mu=mu_new

enddo


!!!!!!!!!!Finish loop

mu=mu_old


write(*,"(a40,f15.10)") "********Mu loop finished, final mu=", mu_new
end subroutine compute_mu_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine compute_nf_kr (nf_kr,nc_kr,b_new_kr,e_kr_tot,w_kr_tot,ef_kr, lambda_new_kr, nfc_kr,diff,nff_kr,mu)
!complex(KIND=idpc) 	:: w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:n_k_kr-1)
!real(KIND=idp) 		:: b_new_kr(0:n_sites_kr-1),lambda_new_kr(0:n_sites_kr-1)
!real(KIND=idp) 		:: E_kr_tot(0:size_h_kr-1,0:n_k_kr-1),diff,mu
!real(kind=idp)		:: nf_kr(0:n_sites_kr-1),nc_kr(0:n_sites_kr-1),nfc_kr(0:n_sites_kr-1),nff_kr(0:n_sites_kr-1)
!real(kind=idp)		:: ec_kr(0:n_sites_kr-1),ef_kr(0:n_sites_kr-1)
!integer		::i,j,ik,is,is2
subroutine compute_nf_r (nf,nc,nfc,nff,b_new,lambda_new,b,lambda,mu,   e_r_tot,w_r_tot,diff, k_vecs_r,nk_r,nk_r_tot, total_size,n_sites,nx,ny,nz,h_psi_tot,index_psi_tot,lvac_r, lambda_par)
integer ,intent(in) 		:: nk_r,nk_r_tot, total_size,n_sites,nx,ny,nz,index_psi_tot(0:total_size-1,n_max_link,3)
complex(KIND=idpc) ,intent(in)	:: w_r_tot(0:total_size-1,0:total_size-1, 0:nk_r-1), h_psi_tot(0:total_size-1,n_max_link,0:nk_r-1)
real(KIND=idp) ,intent(out)	:: b_new(0:n_sites-1),lambda_new(0:n_sites-1)
real(KIND=idp) ,intent(in)	:: b(0:n_sites-1),lambda(0:n_sites-1),mu, lambda_par
real(KIND=idp) ,intent(in) 	:: E_r_tot(0:total_size-1, 0:nk_r-1), k_vecs_r(0:nk_r-1,4)
real(kind=idp) ,intent(out)	:: nf(0:n_sites-1),nc(0:n_sites-1),nfc(0:n_sites-1),nff(0:n_sites-1)
real(kind=idp), allocatable	:: nf_k(:,:),nc_k(:,:),nfc_k(:,:),nff_k(:,:)
real(kind=idp)			:: occup,diff
integer				::i,j,ix,iy,iz,icf,is,o1
integer				::ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8
real(kind=idp) ,allocatable	:: nf_corr(:),nc_corr(:),nfc_corr(:),nff_corr(:)
logical, intent(in)		::lvac_r(0:n_sites-1)


 print *, "nf, nc, etc."
 
!allocate(nf_k(0:n_sites-1, 0:nk_r-1))
!allocate(nc_k(0:n_sites-1, 0:nk_r-1))
!allocate(nfc_k(0:n_sites-1, 0:nk_r-1))
!allocate(nff_k(0:n_sites-1, 0:nk_r-1))
 
 
nf=0
nc=0
nfc=0	!V included
nff=0	!tf included

!nf_k=0
!nc_k=0
!nfc_k=0	!V included
!nff_k=0	!tf included


  do ik=0, nk_r-1  ! k points
  print *, "k=", ik
!  print *, "weight", k_vecs_r(ik,4)/real(nkx_r*nky_r)
   do i=0, total_size-1	!eigenvalues
    occup=1/ (exp((e_r_tot(i,ik)-mu)/T)+1)*k_vecs_r(ik,4)/real(nk_r_tot)
    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z
       do icf=0,1	!c/f
        do is=0,1	!up/down 
         do o1=0,l_orb(icf)-1

	ind=index_xyz(ix,iy,iz,nx,ny)					!site
	indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz)
	
!occupation
!!!!! nf
if (icf==1) nf (ind)= nf (ind) + (abs(w_r_tot(indp,i,ik)))**2 * occup
!if (icf==1) nf_k (ind,ik)= nf_k (ind,ik) + (abs(w_r_tot(indp,i,ik)))**2 * occup
!!!!! nc
if (icf==0) nc (ind)= nc (ind) + (abs(w_r_tot(indp,i,ik)))**2 * occup	
!if (icf==0) nc_k (ind,ik)= nc_k (ind,ik) + (abs(w_r_tot(indp,i,ik)))**2 * occup	

!!!nfc,nff

     do j=1,n_max_link
      ind2 =index_psi_tot(indp,j,1)
      ind2p=index_psi_tot(indp,j,2)
      icf2 =index_psi_tot(indp,j,3)
    
  
if (icf==1 .and. icf2==0) nfc(ind)=nfc(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   *occup
if (icf==1 .and. icf2==1) nff(ind)=nff(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		   	 	            *b(ind2) &
		  	  		   *occup
!nff=2*nff_k*b
!if (icf==1 .and. icf2==0) nfc_k(ind,ik)=nfc_k(ind,ik)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
!		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
!		  	  		   *occup
!if (icf==1 .and. icf2==1) nff_k(ind,ik)=nff_k(ind,ik)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
!		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
!		   	 	            *b(ind2) &
!		  	  		   *occup



	   enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo !ik
   
!correct occupations for reduced set of k points

if (lkred_r .and. nx==ny .and. nx>1 .and. lysameasx .and. .not. ldisorder_r) then

allocate(nc_corr(0:n_sites-1))
allocate(nf_corr(0:n_sites-1))
allocate(nfc_corr(0:n_sites-1))
allocate(nff_corr(0:n_sites-1))

    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z

	ind=ix + iy*nx + iz*nx*ny		!ij
	ind_2=(nx-1-ix) + iy*nx + iz*nx*ny	!-ij
	ind_3=ix + (ny-1-iy)*nx + iz*nx*ny		!i-j
	ind_4=(nx-1-ix) + (ny-1-iy)*nx + iz*nx*ny		!-i-j
	ind_5=iy + ix*nx + iz*nx*ny		!ji
	ind_6=(ny-1-iy) + ix*nx + iz*nx*ny		!-ji
	ind_7=iy + (nx-1-ix)*nx + iz*nx*ny		!j-i
	ind_8=(ny-1-iy) + (nx-1-ix)*nx + iz*nx*ny		!-j-i


nf_corr(ind)=(nf(ind)+nf(ind_2)+nf(ind_3)+nf(ind_4)+nf(ind_5)+nf(ind_6)+nf(ind_7)+nf(ind_8))/8
nc_corr(ind)=(nc(ind)+nc(ind_2)+nc(ind_3)+nc(ind_4)+nc(ind_5)+nc(ind_6)+nc(ind_7)+nc(ind_8))/8
nfc_corr(ind)=(nfc(ind)+nfc(ind_2)+nfc(ind_3)+nfc(ind_4)+nfc(ind_5)+nfc(ind_6)+nfc(ind_7)+nfc(ind_8))/8
nff_corr(ind)=(nff(ind)+nff(ind_2)+nff(ind_3)+nff(ind_4)+nff(ind_5)+nff(ind_6)+nff(ind_7)+nff(ind_8))/8

enddo
enddo
enddo

nf=nf_corr
nc=nc_corr
nfc=nfc_corr
nff=nff_corr

deallocate(nc_corr)
deallocate(nf_corr)
deallocate(nfc_corr)
deallocate(nff_corr)


endif

Nftot=sum(nf(0:n_sites-1))	!total number of f electrons
Nctot=sum(nc(0:n_sites-1))	!total number of c electrons
Nfctot=Nftot+Nctot		!total number of f+c electrons
fctot=sum(nfc(0:n_sites-1))	!total <Vf^+c + V*c^+f>

!!!!!!!!!total free energy
!da fare

!!!!!!


!!!!!!set new variables

do i=0, n_sites-1
!print *, scf_type
if (scf_type=='sb') then
 b_new(i)=(nfc(i)+nff(i))/(2*lambda(i))
 lambda_new(i)=lambda(i)+lambda_par*(1-b_new(i)**2-nf(i))
! ef(i,:)=ef_site(i,:)-lambda_new(i)-mu_new
else if (scf_type=='ga') then
 b_new(i)=(nfc(i)+nff(i))/(4*lambda(i))*(2.-b(i)**2)**2
 lambda_new(i)=lambda(i)+lambda_par*((1-b_new(i)**2)/(1-(b_new(i)**2)/2)-nf(i))
endif

!b_new_k=nfc_k/2/(-nff_k+2*lambda_k/(2-b_k**2)**2)
!lambda_new_k=lambda_k+lambda_par_k*((1-b_new_k**2)/(1-b_new_k**2/2)-nf_k)



!for a vacancy
if (lvac_r(i)) b_new(i)=0
if (lvac_r(i)) lambda_new(i)=-10


enddo



  diff=0

do i=0, n_sites-1
if ( .not. lvac_r(i)) then
if (scf_type=='sb')  diff=max(diff, (1-b_new(i)**2-nf(i))**2+  (b_new(i)-(nfc(i)+nff(i))/(2*lambda_new(i)))**2)
if (scf_type=='ga')  diff=max(diff, ((1-b_new(i)**2)/(1-b_new(i)**2/2)-nf(i))**2+(lambda_new(i)*dpdr(b_new(i))+nfc(i)+nff(i))**2)
endif
enddo

diff=sqrt(diff)

  write(*, "(a10, i5, f20.15, 2(f20.10))"), "?", i_iter, diff, E_tot, E_tot_l
  print *, ""



!!!!!!!!!!Writes new parameters
!write(*,"(3a4,100a15)") "ix","iy","iz", "nf_k", "nc_k", "nfc_k",  "nff_k"
do ik=0, nk_r-1
!write(*, "(i5, 4f15.10)"), ik, k_vecs_r(ik,1), k_vecs_r(ik,2), k_vecs_r(ik,3), k_vecs_r(ik,4)
do ix=0,nx-1	!x
 do iy=0,ny-1	!y
  do iz=0,nz-1	!z

!ind=ix + iy*nx + iz*nx*ny		
ind=index_xyz(ix,iy,iz,nx,ny)
!write(*,"(3i4,100f15.10)") ix,iy,iz, nf_k(ind,ik), nc_k(ind,ik), nfc_k(ind,ik),  nff_k(ind,ik)
  enddo
 enddo
enddo
enddo

print *, " "


write (*, "(3a4,12(a14))"), "x", "y", "z", "b_old","b_new","nf","bncf","b_**2+nf","nfc","lambda","lambda_new", "nff"

do ix=0,nx-1	!x
 do iy=0,ny-1	!y
  do iz=0,nz-1	!z

!ind=ix + iy*nx + iz*nx*ny		
ind=index_xyz(ix,iy,iz,nx,ny)
write(*,"(3i4,100f14.9)") ix,iy,iz, b(ind), b_new(ind), nf(ind),(nfc(ind)+nff(ind))/(2*lambda_new(ind)), (real(b_new(ind)))**2+nf(ind), nfc(ind),lambda(ind), lambda_new(ind), nff(ind)

  enddo
 enddo
enddo


!write(*,"(a10,100f20.10)") "!energy=", E_tot, E_tot_l


write(*,"(10(a15))"), " "
write(*,"(10(a15))") "filling", "ntot", "nf/n_sites","nc/n_sites","nf","nc"  
write(*,"(100f15.10)") , Nfctot/n_sites, Nfctot,Nftot/n_sites, Nctot/n_sites, nftot, nctot
!write(*,"(a50,i5,100f10.5)") "iter_mu, mu_guess, mu=", i_iter_mu, mu_guess, mu
print *, diff


!deallocate(nf_k)
!deallocate(nc_k)
!deallocate(nfc_k)
!deallocate(nff_k)



end subroutine compute_nf_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

complex function sph_arm(l,m,theta,phi)
integer::m,l
real(KIND=idp):: theta,phi

if (l==3) then
if(m==-3) sph_arm=+1.0d0/8	*sqrt(35.0d0/pi)	*exp(-3.0d0*(0,1)*phi)	*sin(theta)**3
if(m==-2) sph_arm=+1.0d0/4	*sqrt(105.0d0/2.0d0/pi)	*exp(-2.0d0*(0,1)*phi)	*sin(theta)**2*cos(theta)
if(m==-1) sph_arm=+1.0d0/8	*sqrt(21.0d0/pi)	*exp(-1.0d0*(0,1)*phi)	*sin(theta)*(5.0d0*cos(theta)**2-1.0d0)
if(m==0)  sph_arm=+1.0d0/4	*sqrt(7.0d0/pi)		*(5.0d0*cos(theta)**3-3.0d0*cos(theta))
if(m==1)  sph_arm=-1.0d0/8	*sqrt(21.0d0/pi)	*exp(1.0d0*(0,1)*phi)	*sin(theta)*(5.0d0*cos(theta)**2-1.0d0)
if(m==2)  sph_arm=+1.0d0/4	*sqrt(105.0d0/2.0d0/pi)	*exp(2.0d0*(0,1)*phi)	*sin(theta)**2*cos(theta)
if(m==3)  sph_arm=-1.0d0/8	*sqrt(35.0d0/pi)	*exp(3.0d0*(0,1)*phi)	*sin(theta)**3

elseif (l==1) then
if(m==-1) sph_arm=+1.0d0/2	*sqrt(3.0d0/2.0d0/pi)	*exp(-1.0d0*(0,1)*phi)	*sin(theta)
if(m==0)  sph_arm=+1.0d0/2	*sqrt(3.0d0/pi)		*cos(theta)
if(m==1)  sph_arm=-1.0d0/2	*sqrt(3.0d0/2.0d0/pi)	*exp(+1.0d0*(0,1)*phi)	*sin(theta)



endif

!sph_arm=ylm(3,m,theta,phi)

return

end function sph_arm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!real function fac10 (n)
!integer ::n
!fac10=1
!end function fac10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!d_z^2 and d_x^2-y^2
complex function sph_arm_d(m,theta,phi)
integer::m
real(KIND=idp):: theta,phi

if(m==1)  sph_arm_d=+1./4*sqrt(5/pi)				*(3*cos(theta)**2-1)		!z^2
if(m==0)  sph_arm_d=+1./4*sqrt(15/pi)	*(cos(phi)**2-sin(phi)**2)*sin(theta)**2		!x^2-y^2
if(m==2)  sph_arm_d=0

return

end function sph_arm_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp)      function CG2(j1,j2,j3,m1,m2,m3)
   DOUBLE PRECISION :: J1, J2, J3, M1, M2, M3

!j1=3
!j2=1-2
!j3=5-2
 
 cg2=0

 if (abs(m1+m2-m3)<1e-6) cg2=2*m2*sqrt((7.0d0/2.0d0-m3*2*m2)/7.0d0)



end function cg2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=idp)      function CG(j1,j2,j3,m1,m2,m3)

      IMPLICIT NONE

      INTEGER :: I, K
      DOUBLE PRECISION :: J1, J2, J3, M1, M2, M3, C, SUMK, TERM
      DOUBLE PRECISION, DIMENSION(0:99) :: FACT

      LOGICAL :: ISFRAC




!
!     Compute table of factorials.
!

      FACT(0) = 1.0D0

      DO I = 1, 99
         FACT(I) = I * FACT(I-1)
      END DO

!
!     Check for invalid input.
!

!      IF (ISFRAC(J1+J2+J3) .OR. ISFRAC(J1+M1)     .OR. ISFRAC(J2+M2) .OR.  &
!          ISFRAC(J3+M3)    .OR. ISFRAC(-J1+J3-M2) .OR. ISFRAC(-J2+J3+M1)) THEN
!         WRITE (UNIT=*, FMT='(/A)') ' Invalid input.'
!         STOP
 !     END IF

!
!     Check for conditions that give C = 0.
!

      IF ( (J3 .LT. ABS(J1-J2)) .OR.  &
           (J3 .GT. (J1+J2))    .OR.  &
           (ABS(M1) .GT. J1)    .OR.  &
           (ABS(M2) .GT. J2)    .OR.  &
           (ABS(M3) .GT. J3)) THEN
         C = 0.0D0
      ELSE

!
!     Compute Clebsch-Gordan coefficient.
!

         C = SQRT((J3+J3+1)/FACT(NINT(J1+J2+J3+1)))
         C = C * SQRT(FACT(NINT(J1+J2-J3))*FACT(NINT(J2+J3-J1))*FACT(NINT(J3+J1-J2)))
         C = C * SQRT(FACT(NINT(J1+M1))*FACT(NINT(J1-M1))*FACT(NINT(J2+M2))*FACT(NINT(J2-M2))*FACT(NINT(J3+M3))*FACT(NINT(J3-M3)))
         SUMK = 0.0D0
         DO K = 0, 99
            IF (J1+J2-J3-K .LT. 0.0D0) CYCLE
            IF (J3-J1-M2+K .LT. 0.0D0) CYCLE
            IF (J3-J2+M1+K .LT. 0.0D0) CYCLE
            IF (J1-M1-K    .LT. 0.0D0) CYCLE
            IF (J2+M2-K    .LT. 0.0D0) CYCLE
            TERM = FACT(NINT(J1+J2-J3-K))*FACT(NINT(J3-J1-M2+K))*FACT(NINT(J3-J2+M1+K))*FACT(NINT(J1-M1-K))*  &
               FACT(NINT(J2+M2-K))*FACT(K)
            IF (MOD(K,2) .EQ. 1) TERM = -TERM
            SUMK = SUMK + 1.0D0/TERM
         END DO
         C = C * SUMK
      END IF

!
!     Print result.
!

!      WRITE (UNIT=*, FMT='(/1X,F15.10)') C

if (abs(m3- m1-m2)>1d-7) c=0.0d0
!write (*,"(10(f10.5))"), m1,m2,m3,abs(m3- m1-m2),c

!      STOP

 cg=c
return
      END function CG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (kind=idp) function lorentz(x,g)
real (kind=idp),intent(in)	::x,g

lorentz=g/pi/(x**2+g**2)
!lorentz=1/g/(exp(x/g)+1)/(exp(-x/g)+1)
return


end function lorentz

!!!

real (kind=idp) function gauss(x,g)
real (kind=idp),intent(in)        ::x,g

gauss=exp(-(x**2/2/g**2))/sqrt(pi)
return

end function gauss

!!!

real (kind=idp) function fermi(x,g)
real (kind=idp),intent(in)        ::x,g

fermi=1/g/(exp(x/g)+1)/(exp(-x/g)+1)
return


end function fermi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_nn(n_nn, n_nn2, n_nn3)
integer,intent(out)	::n_nn, n_nn2, n_nn3
integer 		::n_nn_tot

!simple cubic
n_nn=0
if (lkx) n_nn=n_nn+2
if (lky) n_nn=n_nn+2
if (lkz) n_nn=n_nn+2

n_nn2=0
if (lkx2 .and. lky2) n_nn2=n_nn2+4
if (lky2 .and. lkz2) n_nn2=n_nn2+4
if (lkz2 .and. lkx2) n_nn2=n_nn2+4


n_nn3=0
if (lkx .and. lky .and. lkz) n_nn3=8


!not read
if (.not. lreadham) then
n_max_link=(n_nn+n_nn2+n_nn3)*dim_hk+2	!+2: energy and onsite hybr
n_link_tot=(1+n_nn+n_nn2+n_nn3)*dim_hk
n_nn_tot=n_nn+n_nn2+n_nn3
!read
else
n_link_tot=(2*n_shells+1)**3*dim_hk

print *, "Nns:"
n_nn_tot=0
do lx=-n_shells, n_shells
do ly=-n_shells, n_shells
do lz=-n_shells, n_shells
if (sqrt(real(lx**2+ly**2+lz**2)) <=dist_max .and. sqrt(real(lx**2+ly**2+lz**2))>0.001) then
write (*, "(10i5)"), n_nn_tot+1, lx,ly,lz
n_nn_tot=n_nn_tot+1
endif
enddo
enddo
enddo

n_max_link=n_nn_tot*dim_hk+1	!1: onsite energy

endif


print *, "n_nn:",n_nn
print *, "n_nn2:",n_nn2
print *, "n_nn3:",n_nn3
print *, "n_nn_tot:",n_nn_tot
print *, "n_link_tot:",n_link_tot
print *, "n_max_link:",n_max_link


allocate(index_psi(n_max_link,8))
allocate(h_psi(n_max_link))


end subroutine compute_nn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine set_nn(nn_theta_phi, nn_xyz, nn_ijk)
real(kind=idp)		::nn_theta_phi(n_nn,2), nn_xyz(n_nn,3)
integer			::nn_ijk(n_nn,3)
integer		::i,j

!!!!!!!!!!!!
!correction in tight-binding nn for different surfaces
!surface momenta
!rotation to get the spin in the x,y,z system of the surface
k_surf=0
x_surf=0
nn_ss=0
r1r2r3=0

if (surf_type=="001") then
r1a=1
r1b=0
r1c=0
r2a=0
r2b=1
r2c=0
r3a=0
r3b=0
r3c=1
elseif (surf_type=="100") then
r1a=0
r1b=0
r1c=1
r2a=1
r2b=0
r2c=0
r3a=0
r3b=1
r3c=0
elseif (surf_type=="110") then
r1a=0
r1b=0
r1c=1
r2a=1
r2b=-1
r2c=0
r3a=1
r3b=0
r3c=0
elseif (surf_type=="111") then
r1a=1
r1b=0
r1c=-1
r2a=0
r2b=1
r2c=-1
r3a=0
r3b=0
r3c=1
elseif (surf_type=="210") then
r1a=0
r1b=0
r1c=1
r2a=2
r2b=-1
r2c=0
r3a=1
r3b=0
r3c=0
else	!001
endif


!
! r1 r2 r3 as a function of a=x b=y c=z
r1r2r3(1,1)=r1a
r1r2r3(1,2)=r1b
r1r2r3(1,3)=r1c
r1r2r3(2,1)=r2a
r1r2r3(2,2)=r2b
r1r2r3(2,3)=r2c
r1r2r3(3,1)=r3a
r1r2r3(3,2)=r3b
r1r2r3(3,3)=r3c

det_r1r2r3=det3_int(r1r2r3)

!to correct if the 3rd direction is wrong
if (det_r1r2r3==-1) then
r1r2r3(3,1)=-r3a
r1r2r3(3,2)=-r3b
r1r2r3(3,3)=-r3c
endif

det_r1r2r3=det3_int(r1r2r3)

!plane_index: direction perpendicular to plane defined by r1 and r2 -> vector product
! this serves both as a third direction and as plane index
planeindex(1)=r1r2r3(1,2)*r1r2r3(2,3)-r1r2r3(1,3)*r1r2r3(2,2)
planeindex(2)=r1r2r3(1,3)*r1r2r3(2,1)-r1r2r3(1,1)*r1r2r3(2,3)
planeindex(3)=r1r2r3(1,1)*r1r2r3(2,2)-r1r2r3(1,2)*r1r2r3(2,1)

r1dotr2=r1r2r3(1,1)*r1r2r3(2,1)+r1r2r3(1,2)*r1r2r3(2,2)+r1r2r3(1,3)*r1r2r3(2,3)
r1square=r1r2r3(1,1)**2+r1r2r3(1,2)**2+r1r2r3(1,3)**2
r2square=r1r2r3(2,1)**2+r1r2r3(2,2)**2+r1r2r3(2,3)**2
cosr1r2=r1dotr2/sqrt(dble(r1square))/sqrt(dble(r2square))

!rotation matrix from one to the other ref system
! i normalize the first vector
rot_mat_surf(1,1)=dble(r1r2r3(1,1))/sqrt(dble(r1r2r3(1,1)**2+r1r2r3(1,2)**2+r1r2r3(1,3)**2))
rot_mat_surf(1,2)=dble(r1r2r3(1,2))/sqrt(dble(r1r2r3(1,1)**2+r1r2r3(1,2)**2+r1r2r3(1,3)**2))
rot_mat_surf(1,3)=dble(r1r2r3(1,3))/sqrt(dble(r1r2r3(1,1)**2+r1r2r3(1,2)**2+r1r2r3(1,3)**2))
!i ortogonalize the second vector wrt to the first one
r2_temp(1)=dble(r1r2r3(2,1))-cosr1r2*dble(r1r2r3(1,1))
r2_temp(2)=dble(r1r2r3(2,2))-cosr1r2*dble(r1r2r3(1,2))
r2_temp(3)=dble(r1r2r3(2,3))-cosr1r2*dble(r1r2r3(1,3))
!then i normalize
rot_mat_surf(2,1)=r2_temp(1)/sqrt(r2_temp(1)**2+r2_temp(2)**2+r2_temp(3)**2)
rot_mat_surf(2,2)=r2_temp(2)/sqrt(r2_temp(1)**2+r2_temp(2)**2+r2_temp(3)**2)
rot_mat_surf(2,3)=r2_temp(3)/sqrt(r2_temp(1)**2+r2_temp(2)**2+r2_temp(3)**2)
!the third one is already perpendicular to the first two
rot_mat_surf(3,1)=dble(planeindex(1))/sqrt(dble(planeindex(1)**2+planeindex(2)**2+planeindex(3)**2))
rot_mat_surf(3,2)=dble(planeindex(2))/sqrt(dble(planeindex(1)**2+planeindex(2)**2+planeindex(3)**2))
rot_mat_surf(3,3)=dble(planeindex(3))/sqrt(dble(planeindex(1)**2+planeindex(2)**2+planeindex(3)**2))

det_rotmat=det3_real(rot_mat_surf)

beta_rot=acos(rot_mat_surf(3,3))
!if(rot_mat_surf(1,3)*rot_mat_surf(2,3)<=0) alfa_rot=asin(rot_mat_surf(2,3)/sin(beta_rot))
!if(rot_mat_surf(1,3)*rot_mat_surf(2,3)>0)  alfa_rot=-asin(rot_mat_surf(2,3)/sin(beta_rot))
!if(rot_mat_surf(3,1)*rot_mat_surf(3,2)>=0) gamma_rot=asin(rot_mat_surf(3,2)/sin(beta_rot))
!if(rot_mat_surf(3,1)*rot_mat_surf(3,2)<0)  gamma_rot=-asin(rot_mat_surf(3,2)/sin(beta_rot))
alfa_rot=atan2(rot_mat_surf(2,3)/sin(beta_rot),-rot_mat_surf(1,3)/sin(beta_rot))
gamma_rot=atan2(rot_mat_surf(3,2)/sin(beta_rot),rot_mat_surf(3,1)/sin(beta_rot))

if (surf_type=="111") alfa_rot=-pi/2

!invert r1r2r3 vs abc
r1r2r3_inv(1,1)=r1r2r3(2,2)*r1r2r3(3,3)-r1r2r3(2,3)*r1r2r3(3,2)
r1r2r3_inv(1,2)=r1r2r3(1,3)*r1r2r3(3,2)-r1r2r3(3,3)*r1r2r3(1,2)
r1r2r3_inv(1,3)=r1r2r3(1,2)*r1r2r3(2,3)-r1r2r3(2,2)*r1r2r3(1,3)
r1r2r3_inv(2,1)=r1r2r3(2,3)*r1r2r3(3,1)-r1r2r3(2,1)*r1r2r3(3,3)
r1r2r3_inv(2,2)=r1r2r3(1,1)*r1r2r3(3,3)-r1r2r3(1,3)*r1r2r3(3,1)
r1r2r3_inv(2,3)=r1r2r3(1,3)*r1r2r3(2,1)-r1r2r3(1,1)*r1r2r3(2,3)
r1r2r3_inv(3,1)=r1r2r3(2,1)*r1r2r3(3,2)-r1r2r3(2,2)*r1r2r3(3,1)
r1r2r3_inv(3,2)=r1r2r3(1,2)*r1r2r3(3,1)-r1r2r3(1,1)*r1r2r3(3,2)
r1r2r3_inv(3,3)=r1r2r3(1,1)*r1r2r3(2,2)-r1r2r3(1,2)*r1r2r3(2,1)

!matrix for hoppings
do i=1,3
do j=1,3
nn_ss(i,j)=r1r2r3_inv(j,i)
enddo
enddo


!if r1 and r2 perpendicular: rectangular WSC and BZ
if(r1dotr2==0) then
x_surf(1,1)=sqrt(dble(r1r2r3(1,1)**2+r1r2r3(1,2)**2+r1r2r3(1,3)**2))
x_surf(2,2)=sqrt(dble(r1r2r3(2,1)**2+r1r2r3(2,2)**2+r1r2r3(2,3)**2))
!r1 and r2 same magnitude: (x,y), (-x,y)
elseif(r1square==r2square) then
x_surf(1,1)=sqrt(dble(r1square))*sin(acos(cosr1r2)/2)
x_surf(1,2)=sqrt(dble(r1square))*cos(acos(cosr1r2)/2)
x_surf(2,1)=-x_surf(1,1)
x_surf(2,2)=+x_surf(1,2)
!different: one along x, the other consequently
else
x_surf(1,1)=sqrt(dble(r1square))
x_surf(2,1)=sqrt(dble(r2square))*cosr1r2
x_surf(2,2)=sqrt(dble(r2square))*sqrt(1-cosr1r2**2)
endif


!k_vecs 2d in cartesian coordinates
k_surf(1,1)=x_surf(2,2)/(x_surf(1,1)*x_surf(2,2)-x_surf(1,2)*x_surf(2,1))	!k1x
k_surf(1,2)=-x_surf(2,1)/(x_surf(1,1)*x_surf(2,2)-x_surf(1,2)*x_surf(2,1))	!k1y
k_surf(2,1)=x_surf(1,2)/(x_surf(2,1)*x_surf(1,2)-x_surf(2,2)*x_surf(1,1))	!k1x
k_surf(2,2)=-x_surf(1,1)/(x_surf(2,1)*x_surf(1,2)-x_surf(2,2)*x_surf(1,1))	!k1y




alfa_rot_mat=0
beta_rot_mat=0
gamma_rot_mat=0
rot_mat_temp=0
rot_mat_surf_check=0
!z
alfa_rot_mat(1,1)=cos(alfa_rot)
alfa_rot_mat(1,2)=sin(alfa_rot)
alfa_rot_mat(2,1)=-sin(alfa_rot)
alfa_rot_mat(2,2)=cos(alfa_rot)
alfa_rot_mat(3,3)=1
!y
beta_rot_mat(1,1)=cos(beta_rot)
beta_rot_mat(1,3)=-sin(beta_rot)
beta_rot_mat(2,2)=1
beta_rot_mat(3,1)=sin(beta_rot)
beta_rot_mat(3,3)=cos(beta_rot)
!z
gamma_rot_mat(1,1)=cos(gamma_rot)
gamma_rot_mat(1,2)=sin(gamma_rot)
gamma_rot_mat(2,1)=-sin(gamma_rot)
gamma_rot_mat(2,2)=cos(gamma_rot)
gamma_rot_mat(3,3)=1

do i=1,3
do j=1,3
do k=1,3
rot_mat_temp(i,j)=rot_mat_temp(i,j)+beta_rot_mat(i,k)*gamma_rot_mat(k,j)
enddo
enddo
enddo

do i=1,3
do j=1,3
do k=1,3
rot_mat_surf_check(i,j)=rot_mat_surf_check(i,j)+alfa_rot_mat(i,k)*rot_mat_temp(k,j)
enddo
enddo
enddo

!!!!!!!print
OPEN(unit=300,file=trim(label)//'/surface',status='unknown')
do k=6,300,294

write(k, "(a30)"), "r1r2r3= M abc"
do i=1,3
write(k, "(3i5)"), (r1r2r3(i,j),j=1,3)
enddo

write (k, "(a15, i5)"),"determinant=", det_r1r2r3

write(k, "(a30)"),  "plane indices"
write(k, "(3i5)"), (planeindex(j),j=1,3)

write(k, "(a30)"),  "rotation matrix"
do i=1,3
write(k, "(3f10.5)"), (rot_mat_surf(i,j), j=1,3)
enddo
write(k, "(a15,f10.5)"), "determinant:", det_rotmat

write(k, "(a30)"),  "alfa,beta,gamma"
write(k, "(3f15.10)"), alfa_rot,beta_rot,gamma_rot
write(k, "(a30)"),  "tan alfa,tan beta,tan gamma"
write(k, "(3e15.6)"), tan(alfa_rot),tan(beta_rot),tan(gamma_rot)
write(k, "(a30)"),  "sin alfa,sin beta,sin gamma"
write(k, "(3f15.10)"), sin(alfa_rot),sin(beta_rot),sin(gamma_rot)
write(k, "(a30)"), "cos alfa,cos beta,cos gamma"
write(k, "(3f15.10)"), cos(alfa_rot),cos(beta_rot),cos(gamma_rot)

write(k, "(a30)"),  "abc= M^-1 r1r2r3"
do i=1,3
write(k, "(3i5)"), (r1r2r3_inv(i,j),j=1,3)
enddo

write(k, "(a30)"),  "hopping matrix"
do i=1,3
write(k, "(3i5)"), (nn_ss(i,j),j=1,3)
enddo

write(k, "(a30)"),  "x vec"
write (k, "(a3, 2f10.5)"),"1", x_surf(1,1), x_surf(1,2)
write (k, "(a3, 2f10.5)"),"2", x_surf(2,1), x_surf(2,2)

write(k, "(a30)"),  "k vec"
write (k, "(a3, 2f10.5)"),"1", k_surf(1,1), k_surf(1,2)
write (k, "(a3, 2f10.5)"),"2", k_surf(2,1), k_surf(2,2)



write(k, "(a30)"), "Rot mat for spin on surface (check)"
do i=1,3
write(k, "(3f10.5)"), (rot_mat_surf_check(i,j), j=1,3)
enddo


enddo	!unitfile k
 close(unit=300)

if (det_r1r2r3 .ne. 1) stop


!!!!!!!!!!!!!!!!!!
!sets nearest neighbors
!simple cubic

nn_theta_phi=0
nn_xyz=0
nn_ijk=0

j=1
if(lkx) then
 nn_ijk(j,1)=1
 nn_ijk(j+1,1)=-1
 nn_theta_phi(j,1)=pi/2
 nn_theta_phi(j,2)=0.0d0
 nn_theta_phi(j+1,1)=pi/2
 nn_theta_phi(j+1,2)=pi
 j=j+2
endif

if(lky) then
 nn_ijk(j,2)=1
 nn_ijk(j+1,2)=-1
 nn_theta_phi(j,1)=pi/2
 nn_theta_phi(j,2)=pi/2
 nn_theta_phi(j+1,1)=pi/2
 nn_theta_phi(j+1,2)=-pi/2
 j=j+2
 endif
 
if(lkz) then
 nn_ijk(j,3)=1
 nn_ijk(j+1,3)=-1
 nn_theta_phi(j,1)=0.0d0
 nn_theta_phi(j,2)=0.0d0
 nn_theta_phi(j+1,1)=pi
 nn_theta_phi(j+1,2)=0.0d0
endif



!coordinates
do i=1,n_nn
nn_xyz(i,1)=sin(nn_theta_phi(i,1))*cos(nn_theta_phi(i,2))	!x
nn_xyz(i,2)=sin(nn_theta_phi(i,1))*sin(nn_theta_phi(i,2))	!y
nn_xyz(i,3)=cos(nn_theta_phi(i,1))				!z
enddo


end subroutine set_nn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_nn2(nn2_theta_phi, nn2_xyz, nn2_ijk)
real(kind=idp)		::nn2_theta_phi(n_nn2,2), nn2_xyz(n_nn2,3)
integer			::nn2_ijk(n_nn2,3)
integer		::i,j


!simple cubic

nn2_theta_phi=0
nn2_xyz=0
nn2_ijk=0

j=1
if(lkx2 .and. lky2 ) then
 !+1 +1 0
 nn2_ijk(j,1)=1
 nn2_ijk(j,2)=1
! nn_theta_phi(j,1)=pi/2
! nn_theta_phi(j,2)=pi/4
 !+1 -1 0
 nn2_ijk(j+1,1)=1
 nn2_ijk(j+1,2)=-1
! nn_theta_phi(j+1,1)=pi/2
! nn_theta_phi(j+1,2)=-pi/4
 !-1 +1 0 
 nn2_ijk(j+2,1)=-1
 nn2_ijk(j+2,2)=1
! nn_theta_phi(j+2,1)=pi/2
! nn_theta_phi(j+2,2)=3*pi/4
 !-1 +1 0 
 nn2_ijk(j+3,1)=-1
 nn2_ijk(j+3,2)=-1
! nn_theta_phi(j+3,1)=pi/2
! nn_theta_phi(j+3,2)=3*pi/4

j=j+4
endif

if(lky2 .and. lkz2) then
 !0 +1 +1
 nn2_ijk(j,2)=1
 nn2_ijk(j,3)=1
 !0 +1 -1
 nn2_ijk(j+1,2)=1
 nn2_ijk(j+1,3)=-1
 !0 -1 +1
 nn2_ijk(j+2,2)=-1
 nn2_ijk(j+2,3)=1
 !0 -1 -1
 nn2_ijk(j+3,2)=-1
 nn2_ijk(j+3,3)=-1
j=j+4
 endif

if(lkz2 .and. lkx2) then
 !+1 0 +1
 nn2_ijk(j,1)=1
 nn2_ijk(j,3)=1
 !+1 0 -1
 nn2_ijk(j+1,1)=1
 nn2_ijk(j+1,3)=-1
 !-1 0 +1
 nn2_ijk(j+2,1)=-1
 nn2_ijk(j+2,3)=1
 !-1  0 -1
 nn2_ijk(j+3,1)=-1
 nn2_ijk(j+3,3)=-1
 endif

 
do j=1, n_nn2
enddo



!coordinates
do i=1,n_nn2
nn2_theta_phi(i,1)=acos(nn2_ijk(i,3)/(sqrt(dble((nn2_ijk(i,1)**2+nn2_ijk(i,2)**2+nn2_ijk(i,3)**2)))))
nn2_theta_phi(i,2)=atan2(dble(nn2_ijk(i,2)),dble(nn2_ijk(i,1)))
nn2_xyz(i,1)=sin(nn2_theta_phi(i,1))*cos(nn2_theta_phi(i,2))	!x
nn2_xyz(i,2)=sin(nn2_theta_phi(i,1))*sin(nn2_theta_phi(i,2))	!y
nn2_xyz(i,3)=cos(nn2_theta_phi(i,1))				!z
enddo


end subroutine set_nn2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_nn3(nn3_theta_phi, nn3_xyz, nn3_ijk)
real(kind=idp)		::nn3_theta_phi(n_nn3,2), nn3_xyz(n_nn3,3)
integer			::nn3_ijk(n_nn3,3)
integer		::i,j


!simple cubic

nn3_theta_phi=0
nn3_xyz=0
nn3_ijk=0

j=1
if(lkx .and. lky .and. lkz) then
 !+1 +1 +1
 nn3_ijk(1,1)=1
 nn3_ijk(1,2)=1
 nn3_ijk(1,3)=1
 !1 1 -1
 nn3_ijk(2,1)=1
 nn3_ijk(2,2)=1
 nn3_ijk(2,3)=-1
 !+1 -1 +1
 nn3_ijk(3,1)=1
 nn3_ijk(3,2)=-1
 nn3_ijk(3,3)=1
 !+1 -1 -1
 nn3_ijk(4,1)=1
 nn3_ijk(4,2)=-1
 nn3_ijk(4,3)=-1
 !-1 +1 +1
 nn3_ijk(5,1)=-1
 nn3_ijk(5,2)=1
 nn3_ijk(5,3)=1
 !-1 1 -1
 nn3_ijk(6,1)=-1
 nn3_ijk(6,2)=1
 nn3_ijk(6,3)=-1
 !-1 -1 +1
 nn3_ijk(7,1)=-1
 nn3_ijk(7,2)=-1
 nn3_ijk(7,3)=1
 !-1 -1 -1
 nn3_ijk(8,1)=-1
 nn3_ijk(8,2)=-1
 nn3_ijk(8,3)=-1


!coordinates
do i=1,n_nn3
nn3_theta_phi(i,1)=acos(nn3_ijk(i,3)/(sqrt(dble((nn3_ijk(i,1)**2+nn3_ijk(i,2)**2+nn3_ijk(i,3)**2)))))
nn3_theta_phi(i,2)=atan2(dble(nn3_ijk(i,2)),dble(nn3_ijk(i,1)))
nn3_xyz(i,1)=sin(nn3_theta_phi(i,1))*cos(nn3_theta_phi(i,2))	!x
nn3_xyz(i,2)=sin(nn3_theta_phi(i,1))*sin(nn3_theta_phi(i,2))	!y
nn3_xyz(i,3)=cos(nn3_theta_phi(i,1))				!z
enddo

endif

end subroutine set_nn3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_index_k(ind, ikx,iky,ikz)
integer,intent(out)	::ikx,iky,ikz
integer,intent(in)	::ind


ikx=mod(ind,nkx_k)
iky=mod((ind-ikx)/nkx_k,nky_k)
ikz=mod((ind-ikx-iky*nkx_k)/nkx_k/nky_k,nkz)

!ind=ix + iy*nx + iz*nx*ny

end subroutine convert_index_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine allocate_var_doublet

allocate(ec_site_k(0:l_orbital-1))
allocate(ef_site_k(0:l_orbital_f-1))
allocate(alpha_jzi(0:l_orbital_f-1,0:1,0:5))
allocate(alpha_jzi_mat(0:5,0:5))
allocate(UtU(0:5,0:5))
allocate(alpha_msigma(0:l_orbital_f-1,0:1,0:1,-3:3))
allocate(phi_sigma_alpha(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1, n_nn))
allocate(phi_alpha_sigma(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1, n_nn))
allocate(phi_sigma_alpha2(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1, n_nn2))
allocate(phi_alpha_sigma2(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1, n_nn2))
allocate(phi_sigma_alpha3(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1, n_nn3))
allocate(phi_alpha_sigma3(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1, n_nn3))
allocate(kin_energy     (0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1, n_nn,0:1))
allocate(kin_energy2     (0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1, n_nn2,0:1))
allocate(kin_energy3     (0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1, n_nn3,0:1))
allocate(nn_ijk(N_nn,3))
allocate(nn_ss(3,3))
allocate(nn_xyz(N_nn,3))
allocate(nn_theta_phi(N_nn,2))
allocate(nn2_ijk(N_nn2,3))
allocate(nn2_xyz(N_nn2,3))
allocate(nn2_theta_phi(N_nn2,2))
allocate(nn3_ijk(N_nn3,3))
allocate(nn3_xyz(N_nn3,3))
allocate(nn3_theta_phi(N_nn3,2))
!allocate(k_vecs(0:n_k-1,1:3))

end subroutine allocate_var_doublet

subroutine deallocate_var_doublet

deallocate(alpha_jzi)
deallocate(UtU)
deallocate(alpha_msigma)
deallocate(phi_sigma_alpha)
deallocate(phi_alpha_sigma)
deallocate(phi_sigma_alpha2)
deallocate(phi_alpha_sigma2)
deallocate(phi_sigma_alpha3)
deallocate(phi_alpha_sigma3)
deallocate(kin_energy)
deallocate(kin_energy2)
deallocate(kin_energy3)
deallocate(nn_ijk)
deallocate(nn_xyz)
deallocate(nn_theta_phi)
deallocate(nn2_ijk)
deallocate(nn2_xyz)
deallocate(nn2_theta_phi)
deallocate(nn3_ijk)
deallocate(nn3_xyz)
deallocate(nn3_theta_phi)
!deallocate(k_vecs)

end subroutine deallocate_var_doublet



subroutine allocate_var_k

!print *, nk3d_irr, "all"
allocate(phi_k_sigma_alpha_tot(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk3d-1))
allocate(phi_k_alpha_sigma_tot(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk3d-1))
allocate(kin_energy_k_tot     (0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk3d-1,0:1))
!allocate(phi_k_sigma_alpha(0:1,0:1))
!allocate(phi_k_alpha_sigma(0:1,0:1))
allocate(ham_k(0:dim_hk-1,0:dim_hk-1))
allocate(e_k(0:dim_hk-1))
allocate(e_k_tot(0:dim_hk-1,0:nk3d-1))
allocate(w_k(0:dim_hk-1,0:dim_hk-1))
allocate(w_k_tot(0:dim_hk-1,0:dim_hk-1,0:nk3d-1))
allocate(V_k_tot(0:dim_hk-1,0:dim_hk-1,0:nk3d-1))
allocate(V_k(0:dim_hk-1,0:dim_hk-1))
!allocate(factor_k_tot(0:nk3d-1))
allocate(delta_k_tot(0:nk3d-1))


end subroutine allocate_var_k

subroutine deallocate_var_k

deallocate(phi_k_sigma_alpha_tot)
deallocate(phi_k_alpha_sigma_tot)
!deallocate(phi_k_sigma_alpha)
!deallocate(phi_k_alpha_sigma)
deallocate(ham_k)
deallocate(e_k)
deallocate(e_k_tot)
deallocate(w_k)
deallocate(w_k_tot)
deallocate(v_k_tot)
deallocate(V_k)
deallocate(k_vecs3d)
end subroutine deallocate_var_k



subroutine allocate_var_r

!allocate(ec(0:n_sites-1))
allocate(ef(0:n_sites-1,0:1))
allocate(ef_site(0:n_sites-1,0:l_orbital_f-1))
allocate(ec_site(0:n_sites-1,0:l_orbital-1))
!allocate(tc(0:n_sites-1,0:n_sites-1))
!allocate(tf(0:n_sites-1,0:n_sites-1))
!allocate(Vcf(0:4*n_sites-1,0:4*n_sites-1))	!sigma and alpha dependent!
allocate(b(0:n_sites-1))
allocate(one(0:n_sites-1))
allocate(zero(0:n_sites-1))
allocate(b_new(0:n_sites-1))
allocate(lambda(0:n_sites-1))
allocate(lambda_new(0:n_sites-1))
allocate(ham(0:total_size-1,0:total_size-1))
allocate(e(0:total_size-1))
allocate(w(0:total_size-1,0:total_size-1))
allocate(nf(0:n_sites-1))	
allocate(nc(0:n_sites-1))
allocate(nfc(0:n_sites-1))
allocate(nff(0:n_sites-1))
if (n_sites <5) allocate(ham_rot(0:total_size-1,0:total_size-1))
if (n_sites <5) allocate(U(0:total_size-1,0:total_size-1))
if (n_sites <5) allocate(Um1(0:total_size-1,0:total_size-1))
if (n_sites <5) allocate(UstarU(0:total_size-1,0:total_size-1))
allocate(index_psi_tot(0:total_size-1,n_max_link,3))
allocate(lvac_r(0:n_sites-1))


end subroutine allocate_var_r

subroutine deallocate_var_r

!allocate(ec(0:n_sites-1))
deallocate(ef)
deallocate(ef_site)
deallocate(ec_site)
!allocate(tc(0:n_sites-1,0:n_sites-1))
!allocate(tf(0:n_sites-1,0:n_sites-1))
!allocate(Vcf(0:4*n_sites-1,0:4*n_sites-1))	!sigma and alpha dependent!
deallocate(b)
deallocate(one)
deallocate(zero)
deallocate(b_new)
deallocate(lambda)
deallocate(lambda_new)
deallocate(ham)
deallocate(e)
deallocate(w)
deallocate(nf)	
deallocate(nc)
deallocate(nfc)
deallocate(nff)
if (n_sites <5) deallocate(ham_rot)
if (n_sites <5) deallocate(U)
if (n_sites <5) deallocate(Um1)
if (n_sites <5) deallocate(UstarU)
!deallocate(index_psi)
!if(allocated(h_psi)) deallocate(h_psi)
if(allocated(index_psi_tot))deallocate(index_psi_tot)
if(allocated(h_psi_tot))deallocate(h_psi_tot)
deallocate(lvac_r)

if (allocated(w_r_tot)) deallocate(w_r_tot)
if (allocated(e_r_tot))deallocate(e_r_tot)

end subroutine deallocate_var_r



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine open_files

!if (loopr .and. ldos) OPEN(unit=10,file=trim(label)//'/dos_r',status='unknown')
if (loopr) OPEN(unit=11,file=trim(label)//'/b_r',status='unknown')
if (loopr) OPEN(unit=51,file=trim(label)//'/b_r0',status='unknown')
!if (loopk .and. ldos) OPEN(unit=12,file=trim(label)//'/dos_k',status='unknown')
if (loopk) OPEN(unit=13,file=trim(label)//'/b_k',status='unknown', position='append')
if (loopk) OPEN(unit=14,file=trim(label)//'/bands_k',status='unknown')
!if (loopkr .and. ldos) OPEN(unit=20,file=trim(label)//'/dos_kr',status='unknown')
if (loopkr) OPEN(unit=21,file=trim(label)//'/bands_kr',status='unknown')
if (loopkr) OPEN(unit=22,file=trim(label)//'/b_kr',status='unknown')
!if (loopkr .and. lwfc) OPEN(unit=23,file=trim(label)//'/wfc_kr',status='unknown')
!if(.not. lparamopen) OPEN(unit=24,file=trim(label)//'/param_info',status='unknown')
!if (loopk) OPEN(unit=25,file=trim(label)//'/bands_k_interp',status='unknown')
if (loopk .and. lnscfk .and. lpathnscf) OPEN(unit=26,file=trim(label)//'/bands_k_nscf',status='unknown')
!if (loopk .and. lnscfk .and. ldos) OPEN(unit=27,file=trim(label)//'/dos_k_nscf',status='unknown')
!if (loopkr .and. lnscfk .and. lpathnscf) OPEN(unit=28,file=trim(label)//'/bands_kr_nscf',status='unknown')
!if (loopkr .and. lnscfk .and. ldos) OPEN(unit=29,file=trim(label)//'/dos_kr_nscf',status='unknown')
!if (loopk) OPEN(unit=30,file=trim(label)//'/b_k2',status='unknown')
!if (loopkr) OPEN(unit=31,file=trim(label)//'/b_kr2',status='unknown')
!OPEN(unit=32,file=trim(label)//'/param',status='unknown')
if (loopk .and. lnscfk.and. .not. lpathnscf) OPEN(unit=33,file=trim(label)//'/ef_k_nscf',status='unknown')
if (loopkr .and. lnscfk.and. .not. lpathnscf) OPEN(unit=34,file=trim(label)//'/ef_kr_nscf',status='unknown')
if (loopk .and. lnscfk .and. .not. lpathnscf) OPEN(unit=35,file=trim(label)//'/bands_k_nscfa',status='unknown')
if (loopkr .and. lnscfk .and. .not. lpathnscf) OPEN(unit=36,file=trim(label)//'/bands_kr_nscfa',status='unknown')
!if (loopr) OPEN(unit=37,file=trim(label)//'/b_r2',status='unknown')
!38 =wfcr
end subroutine open_files 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine close_files

if (loopr) close(unit=10)
if (loopr) close(unit=11)
if (loopk) close(unit=12)
if (loopk) close(unit=13)
if (loopk) close(unit=14)
if (loopr) close(unit=15)
if (loopr) close(unit=16)
if (loopr) close(unit=17)
if (loopr) close(unit=18)
if (loopr) close(unit=19)
if (loopkr) close(unit=20)
if (loopkr) close(unit=21)
if (loopkr) close(unit=22)
if (loopkr) close(unit=23)
 close(unit=24)
!if (loopk) close(unit=25)
if (loopk .and. lnscfk) close(unit=26)
if (loopk .and. lnscfk) close(unit=27)
if (loopkr .and. lnscfk) close(unit=28)
if (loopkr .and. lnscfk) close(unit=29)

end subroutine close_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_phi !(phi_sigma_alpha, phi_alpha_sigma)
!complex(KIND=idpc):: phi_sigma_alpha(0:l_orbital-1,0:l_orbital-1,0:1,0:1,n_nn),phi_alpha_sigma(0:l_orbital-1,0:l_orbital-1,0:1,0:1,n_nn)
integer		::ikx,iky,ikz


phi_sigma_alpha=0
phi_alpha_sigma=0



  do i=1, n_nn
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0, l_orbital-1		!d
      do o2=0, l_orbital-1		!f for now only gamma8
    
if (lphioverlap) then
      do k=-3,3
       phi_sigma_alpha(0,0,is,is2,i)=phi_sigma_alpha(0,0,is,is2,i)+alpha_msigma(0,is2,is,k)*sph_arm(j1,k,nn_theta_phi(i,1),nn_theta_phi(i,2))*sqrt(pi/3)*eta_v
       phi_alpha_sigma(0,0,is,is2,i)=phi_alpha_sigma(0,0,is,is2,i)+alpha_msigma(0,is,is2,k)*conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*sqrt(pi/3)*eta_v
      enddo


     elseif(phi_type=="smb" .or. phi_type=="sm7") then
if (.not. l_k_dz) then
      do k=-3,3	
      									!	     alpha s m
       phi_sigma_alpha(o1,o2,is,is2,i)=phi_sigma_alpha(o1,o2,is,is2,i)+alpha_msigma(o2,is2,is,k)*sph_arm(j1,k,nn_theta_phi(i,1),nn_theta_phi(i,2))*&
       				sph_arm_d(o1,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2))*4*pi/sqrt(15.0d0)*eta_vz
       phi_alpha_sigma(o1,o2,is,is2,i)=phi_alpha_sigma(o1,o2,is,is2,i)+alpha_msigma(o1,is,is2,k)*&
       				conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*sph_arm_d(o2,nn_theta_phi(i,1),nn_theta_phi(i,2))*4*pi/sqrt(15.0d0)*eta_vz
      enddo
 
 else
 !x
if(i<3 .and. is.ne.is2 .and. o1==0 .and. o2==0)      phi_sigma_alpha(o1,o2,is,is2,i)=+nn_ijk(i,1)*(3.0d0*eta_vz+eta_v)/(4.0d0)
if(i<3 .and. is.ne.is2 .and. o1.ne.o2)   	     phi_sigma_alpha(o1,o2,is,is2,i)=-nn_ijk(i,1)*dsqrt(3.0d0)*(eta_vz-eta_v)/(4.0d0)
if(i<3 .and. is.ne.is2 .and. o1==1 .and. o2==1)      phi_sigma_alpha(o1,o2,is,is2,i)=+nn_ijk(i,1)*(eta_vz+3.0d0*eta_v)/(4.0d0)

!if(i<3 .and. is.ne.is2 .and. o1==0 .and. o2==0)      phi_alpha_sigma(o2,o1,is,is2,i)=-nn_ijk(i,1)*(3*eta_vz+eta_v)/(4.0d0)
!if(i<3 .and. is.ne.is2 .and. o1.ne.o2)    	     phi_alpha_sigma(o2,o1,is,is2,i)=+nn_ijk(i,1)*dsqrt(3.0d0)*(eta_vz-eta_v)/(4.0d0)
!if(i<3 .and. is.ne.is2 .and. o1==1 .and. o2==1)      phi_alpha_sigma(o2,o1,is,is2,i)=-nn_ijk(i,1)*(eta_vz+3*eta_v)/(4.0d0)



!y 
if(i>2 .and. i<5 .and. is.ne.is2 .and. o1==0 .and. o2==0)    phi_sigma_alpha(o1,o2,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*(3.0d0*eta_vz+eta_v)/(4.0d0)*(0,1)
if(i>2 .and. i<5 .and. is.ne.is2 .and. o1.ne.o2)  	     phi_sigma_alpha(o1,o2,is,is2,i)=-nn_ijk(i,2)*sqrt(3.0d0)/(4.0d0)*(2*is-1)*(eta_vz-eta_v)*(0,1)
if(i>2 .and. i<5 .and. is.ne.is2 .and. o1==1 .and. o2==1)    phi_sigma_alpha(o1,o2,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*(eta_vz+3.0d0*eta_v)/(4.0d0)*(0,1)
 
!if(i>2 .and. i<5 .and. is.ne.is2 .and. o1==0 .and. o2==0)    phi_alpha_sigma(o2,o1,is,is2,i)=nn_ijk(i,2)*(2*is-1)*(3*eta_vz+eta_v)/(4.0d0)*(0,1)
!if(i>2 .and. i<5 .and. is.ne.is2 .and. o1.ne.o2)  	     phi_alpha_sigma(o2,o1,is,is2,i)=nn_ijk(i,2)*sqrt(3.0d0)/(4.0d0)*(2*is-1)*(eta_vz-eta_v)*(0,1)
!if(i>2 .and. i<5 .and. is.ne.is2 .and. o1==1 .and. o2==1)    phi_alpha_sigma(o2,o1,is,is2,i)=nn_ijk(i,2)*(2*is-1)*(eta_vz+3*eta_v)/(4.0d0)*(0,1)

 
 
!z
if(i>4 .and. is==is2 .and. o1==1 .and. o2==1) phi_sigma_alpha(o1,o2,is,is2,i)=nn_ijk(i,3)*(2*is-1)*eta_vz	!d_z^2, Gamma_8(2)
if(i>4 .and. is==is2 .and. o1==0 .and. o2==0) phi_sigma_alpha(o1,o2,is,is2,i)=nn_ijk(i,3)*(2*is-1)*eta_v	!d_x^2-y^2, Gamma_8(1)
!if(i>4 .and. is==is2 .and. o1==1 .and. o2==1) phi_alpha_sigma(o2,o1,is,is2,i)=-nn_ijk(i,3)*(2*is-1)*eta_vz
!if(i>4 .and. is==is2 .and. o1==0 .and. o2==0) phi_alpha_sigma(o2,o1,is,is2,i)=-nn_ijk(i,3)*(2*is-1)*eta_v

phi_alpha_sigma(o2,o1,is2,is,i)=-conjg(phi_sigma_alpha(o1,o2,is,is2,i))

endif






     
     elseif(phi_type=="dze") then
!x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/2.*(0,1)
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)/2.*(0,1)
!y
if(i>2 .and. i<5 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)*(2*is-1)/2.*(0,1)
if(i>2 .and. i<5 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,2)*(2*is-1)/2.*(0,1)
!z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)*(2*is-1)*(0,1)/2
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,3)*(2*is-1)*(0,1)/2

!-sigma_x sin x - sigma_y  sin_y + 2 sigma_z sin_z
     elseif(phi_type=="sxy") then
!x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1)/4.
!if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1)/4.*(0,1)*(2*is-1)
!if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/4.
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))
!y
if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)/4.*(2*is-1)*(0,1)
!if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.
!if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.*(0,1)*(2*is-1)
if(i>2 .and. i<5 .and. is.ne.is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))
!z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)*(2*is-1)/2.
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))

     elseif(phi_type=="sm2") then	!sigma_x sin x + sigma_y  sin_y + sigma_z sin_z
!x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/4.*eta_v
!if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1)/4.*(0,1)*(2*is-1)
!if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/4.
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))
!y
if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.*(2*is-1)*(0,1)*eta_v
!if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.
!if(i>2 .and. i<5 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.*(0,1)*(2*is-1)
if(i>2 .and. i<5 .and. is.ne.is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))
!z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)*(2*is-1)/4.*eta_v
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is2,is,i)=-conjg(phi_sigma_alpha(0,0,is,is2,i))


     elseif(phi_type=="dz2") then
!x
if(i<3 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/2.*(0,1)*(2*is-1)
if(i<3 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)/2.*(0,1)*(2*is-1)
!y
if(i>2 .and. i<5 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)*(2*is-1)/2.*(0,1)
if(i>2 .and. i<5 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,2)*(2*is-1)/2.*(0,1)
!z
if(i>4 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)*(0,1)/2
if(i>4 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,3)*(0,1)/2

     elseif(phi_type=="ide") then		!prop to identity,odd phi(-k)=-phi(k) (s_x+x_y+s_z)sigma_0
if(is==is2)    phi_sigma_alpha(0,0,is,is2,i)=(nn_ijk(i,1)+nn_ijk(i,2)+nn_ijk(i,3))
if(is==is2)    phi_alpha_sigma(0,0,is,is2,i)=-(nn_ijk(i,1)+nn_ijk(i,2)+nn_ijk(i,3))

     elseif(phi_type=="trv") then		!phi(-k)=phi(k), non spin dependent (always for nn) sigma_0(c_x+c_y+c_z)
if(is==is2)    phi_sigma_alpha(0,0,is,is2,i)=eta_v/4.
if(is==is2)    phi_alpha_sigma(0,0,is,is2,i)=eta_v/4.

     elseif(phi_type=="six") then		!prop to sigma x
if(is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=(nn_ijk(i,1)+nn_ijk(i,2)+nn_ijk(i,3))!*(2*is-1)
if(is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=-(nn_ijk(i,1)+nn_ijk(i,2)+nn_ijk(i,3))!*(2*is-1)

     elseif(phi_type=="kim") then		!PRB 85, 125128
!x: -i sin(kx) sigma_x-->(100)=-sigma_x (-100) sigma_x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1)/4.0d0 !* (4*0.24430)
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)/4.0d0!* (4*0.24430)
!y:  i sin(ky) sigma_y-->(010)=sigma_y /(0-10) -sigma_y
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)/4.0d0*(0,1)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.0d0*(0,1)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,2)/4.0d0*(0,1)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=+nn_ijk(i,2)/4.0d0*(0,1)!* (4*0.24430)
!z  2i sin(kz) sigma_z-->(001)=2*sigma_z (00-1)=-2sigma_z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=+nn_ijk(i,3)*(2*is-1)/2.0d0!* (4*0.24430)
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,3)*(2*is-1)/2.0d0!* (4*0.24430)

     elseif(phi_type=="kyz") then		!PRB 85, 125128, change k_y<-->k_z
!x: -2i sin(kx) sigma_x-->(100)=-sigma_x (-100) sigma_x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1)/4.0d0 !* (4*0.24430)
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)/4.0d0!* (4*0.24430)
!z:  2i sin(kz) sigma_y-->(001)=sigma_y /(00-1) -sigma_y
if(i>4 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,3)/4.0d0*(0,1)!* (4*0.24430)
if(i>4 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)/4.0d0*(0,1)!* (4*0.24430)
if(i>4 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,3)/4.0d0*(0,1)!* (4*0.24430)
if(i>4 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,3)/4.0d0*(0,1)!* (4*0.24430)
!y  2i sin(ky) sigma_z-->(010)=2*sigma_z (0-10)=-2sigma_z
if(i>2 .and. i<5  .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)*(2*is-1)/2.0d0!* (4*0.24430)
if(i>2 .and. i<5  .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,2)*(2*is-1)/2.0d0!* (4*0.24430)

     elseif(phi_type=="ki1") then		!PRB 85, 125128
!x: -i sin(kx) sigma_x-->(100)=-sigma_x (-100) sigma_x     x sqrt(3)
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1)/4.0d0*sqrt(3.0d0) !* (4*0.24430)
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,1)/4.0d0*sqrt(3.0d0)!* (4*0.24430)
!y:  i sin(ky) sigma_y-->(010)=sigma_y /(0-10) -sigma_y    x sqrt(3)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)/4.0d0*(0,1)*sqrt(3.0d0)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)/4.0d0*(0,1)*sqrt(3.0d0)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,2)/4.0d0*(0,1)*sqrt(3.0d0)!* (4*0.24430)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=+nn_ijk(i,2)/4.0d0*(0,1)*sqrt(3.0d0)!* (4*0.24430)
!z  0
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=0!* (4*0.24430)
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=0!* (4*0.24430)

     elseif(phi_type=="sra") then		!PRL 108, 126807
!x
if(i<3 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,1) 
if(i<3 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1) 
if(i<3 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=+nn_ijk(i,1) 
if(i<3 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,1)
!y
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,2)*(0,1.0d0)
!
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,3)
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,3)

     elseif(phi_type=="srb" .or. phi_type=="src") then		
!x: -i sin(kx) sigma_x-->(100)=-sigma_x (-100) sigma_x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1) 
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)
!y:  i sin(ky) sigma_y-->(010)=sigma_y /(0-10) -sigma_y
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=+nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=+nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,2)*(0,1.0d0)
!z  2i sin(kz) sigma_z-->(001)=2*sigma_z (00-1)=-sigma_z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=+nn_ijk(i,3)*(2*is-1)
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,3)*(2*is-1)

!same as srb, with 2*z
     elseif(phi_type=="sbz") then		
!x: -i sin(kx) sigma_x-->(100)=-sigma_x (-100) sigma_x
if(i<3 .and. is.ne.is2)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,1) 
if(i<3 .and. is.ne.is2)    phi_alpha_sigma(0,0,is,is2,i)=nn_ijk(i,1)
!y:  i sin(ky) sigma_y-->(010)=sigma_y /(0-10) -sigma_y
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_sigma_alpha(0,0,is,is2,i)=nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_sigma_alpha(0,0,is,is2,i)=-nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==0 .and. is2==1)    phi_alpha_sigma(0,0,is,is2,i)=-nn_ijk(i,2)*(0,1.0d0)
if(i>2 .and. i<5 .and. is==1 .and. is2==0)    phi_alpha_sigma(0,0,is,is2,i)=+nn_ijk(i,2)*(0,1.0d0)
!z  2i sin(kz) sigma_z-->(001)=2*sigma_z (00-1)=-2*sigma_z
if(i>4 .and. is==is2)    phi_sigma_alpha(0,0,is,is2,i)=+2.0d0*nn_ijk(i,3)*(2*is-1)
if(i>4 .and. is==is2)    phi_alpha_sigma(0,0,is,is2,i)=-2.0d0*nn_ijk(i,3)*(2*is-1)


     endif
    enddo
   enddo
  enddo
 enddo

!+i --> - sin
!+1 --> i sin

!Gamma_7
if(phi_type=="sm7") then
   do is=0,1	!sigma	!alpha
    do is2=0,1
     do o1=0, l_orbital-1		!d

!if ( .not. l7z) then
!with x^2-y^2
!x
if(i<3 .and. is.ne.is2 .and. o1==0)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,1)*eta_v7/2
if(i<3 .and. is.ne.is2 .and. o1==1)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,1)*sqrt(3.)*eta_v7/2
if(i<3 .and. is.ne.is2 .and. o1==0)      	     phi_alpha_sigma(2,o1,is,is2,i)=+nn_ijk(i,1)*eta_v7/2
if(i<3 .and. is.ne.is2 .and. o1==1)      	     phi_alpha_sigma(2,o1,is,is2,i)=+nn_ijk(i,1)*sqrt(3.)*eta_v7/2
!y ! (check overall sign)
if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==0)      	     phi_sigma_alpha(o1,2,is,is2,i)=+nn_ijk(i,2)*(2*is-1)*eta_v7/2*(0,1)
if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==1)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*sqrt(3.)*eta_v7/2*(0,1)
if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==0)      	     phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*eta_v7/2*(0,1)
if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==1)      	     phi_alpha_sigma(2,o1,is,is2,i)=+nn_ijk(i,2)*(2*is-1)*sqrt(3.)*eta_v7/2*(0,1)
!z
if(i>4 .and. is==is2 .and. o1==0) phi_sigma_alpha(o1,2,is,is2,i)=nn_ijk(i,3)*(2*is-1)*eta_v7	
if(i>4 .and. is==is2 .and. o1==0) phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,3)*(2*is-1)*eta_v7

!else
!with z^2
!if(i<3 .and. is.ne.is2 .and. o1==1)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,1)*(0,1)*eta_v7/4
!if(i<3 .and. is.ne.is2 .and. o1==0)      	     phi_sigma_alpha(o1,2,is,is2,i)=+nn_ijk(i,1)*(0,1)*sqrt(3.)*eta_v7/4
!if(i<3 .and. is.ne.is2 .and. o1==1)      	     phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,1)*(0,1)*eta_v7/4
!if(i<3 .and. is.ne.is2 .and. o1==0)      	     phi_alpha_sigma(2,o1,is,is2,i)=+nn_ijk(i,1)*(0,1)*sqrt(3.)*eta_v7/4


!if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==1)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*eta_v7/4
!if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==0)      	     phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*sqrt(3.)*eta_v7/4
!if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==1)      	     phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*eta_v7/4
!if(i>2 .and. i<5  .and. is.ne.is2 .and. o1==0)      	     phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,2)*(2*is-1)*sqrt(3.)*eta_v7/4

 
!if(i>4 .and. is==is2 .and. o1==1) phi_sigma_alpha(o1,2,is,is2,i)=-nn_ijk(i,3)*(2*is-1)*(0,1)*eta_v7	!d_z^2, Gamma_8(2)
!if(i>4 .and. is==is2 .and. o1==1) phi_alpha_sigma(2,o1,is,is2,i)=-nn_ijk(i,3)*(2*is-1)*(0,1)*eta_v7
!endif

   enddo
    enddo
     enddo

endif	!sm7

enddo


!!!!!!!!!!!!!!!!!2nd nn

phi_sigma_alpha2=0
phi_alpha_sigma2=0

  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0, l_orbital-1		!d
      do o2=0, l_orbital-1		!f for now only gamma8
    
     if (lphioverlap) then
      do k=-3,3
       phi_sigma_alpha2(0,0,is,is2,i)=phi_sigma_alpha2(0,0,is,is2,i)+alpha_msigma(0,is2,is,k)*sph_arm(j1,k,nn2_theta_phi(i,1),nn2_theta_phi(i,2))*sqrt(pi/3)*eta_v2
       phi_alpha_sigma2(0,0,is,is2,i)=phi_alpha_sigma2(0,0,is,is2,i)+alpha_msigma(0,is,is2,k)*conjg(sph_arm(j1,k,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2)))*sqrt(pi/3)*eta_v2
      enddo
endif
enddo
enddo
enddo
enddo
enddo


if (l_orbital==2 .and. l_orbital_f>1) then
  do i=1, n_nn2
   do is=0,1	!sigma	!alpha


if (i<5) then!110
phi_sigma_alpha2(0,0,is,1-is,i)=+etav2_x1*nn2_ijk(i,1)
phi_sigma_alpha2(0,1,is,1-is,i)=+etav2_x2*nn2_ijk(i,1)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,1-is,i)=+etav2_x7*nn2_ijk(i,1)
phi_sigma_alpha2(1,0,is,1-is,i)=+etav2_z1*nn2_ijk(i,1)
phi_sigma_alpha2(1,1,is,1-is,i)=+etav2_z2*nn2_ijk(i,1)
if (l_orbital_f>2)  phi_sigma_alpha2(1,2,is,1-is,i)=+etav2_z7*nn2_ijk(i,1)
phi_sigma_alpha2(0,0,is,1-is,i)=phi_sigma_alpha2(0,0,is,1-is,i)-etav2_x1*(0,1)*(2*is-1)*nn2_ijk(i,2)
phi_sigma_alpha2(0,1,is,1-is,i)=phi_sigma_alpha2(0,1,is,1-is,i)+etav2_x2*(0,1)*(2*is-1)*nn2_ijk(i,2)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,1-is,i)=phi_sigma_alpha2(0,2,is,1-is,i)-etav2_x7*(0,1)*(2*is-1)*nn2_ijk(i,2)
phi_sigma_alpha2(1,0,is,1-is,i)=phi_sigma_alpha2(1,0,is,1-is,i)+etav2_z1*(0,1)*(2*is-1)*nn2_ijk(i,2)
phi_sigma_alpha2(1,1,is,1-is,i)=phi_sigma_alpha2(1,1,is,1-is,i)-etav2_z2*(0,1)*(2*is-1)*nn2_ijk(i,2)
if (l_orbital_f>2) phi_sigma_alpha2(1,2,is,1-is,i)=phi_sigma_alpha2(1,2,is,1-is,i)+etav2_z7*(0,1)*(2*is-1)*nn2_ijk(i,2)
endif

if (i>4 .and. i<9) then	!011
phi_sigma_alpha2(0,0,is,is,i)=(+etav2_x1+sqrt(3.)*etav2_x2+sqrt(3.)*etav2_z1+3*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(0,1,is,is,i)=(+sqrt(3.)*etav2_x1-etav2_x2+3*etav2_z1-sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,is,i)=(-etav2_x7-sqrt(3.)*etav2_z7)/2.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(1,0,is,is,i)=(+sqrt(3.)*etav2_x1+3*etav2_x2-etav2_z1-sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(1,1,is,is,i)=(+3*etav2_x1-sqrt(3.)*etav2_x2-sqrt(3.)*etav2_z1+etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
if (l_orbital_f>2) phi_sigma_alpha2(1,2,is,is,i)=(-sqrt(3.)*etav2_x7+etav2_z7)/2.*(2*is-1)*nn2_ijk(i,3)

phi_sigma_alpha2(0,0,is,1-is,i)=-(etav2_x1-sqrt(3.)*etav2_x2-sqrt(3.)*etav2_z1+3*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,2)*(0,1)
phi_sigma_alpha2(0,1,is,1-is,i)=-(sqrt(3.)*etav2_x1+etav2_x2-3*etav2_z1-sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,2)*(0,1)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,1-is,i)=-(-etav2_x7+sqrt(3.)*etav2_z7)/2.*(2*is-1)*nn2_ijk(i,2)*(0,1)
phi_sigma_alpha2(1,0,is,1-is,i)=-(sqrt(3.)*etav2_x1-3*etav2_x2+etav2_z1-sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,2)*(0,1)
phi_sigma_alpha2(1,1,is,1-is,i)=-(3*etav2_x1+sqrt(3.)*etav2_x2+sqrt(3.)*etav2_z1+etav2_z2)/4.*(2*is-1)*nn2_ijk(i,2)*(0,1)
if (l_orbital_f>2) phi_sigma_alpha2(1,2,is,1-is,i)=-(-sqrt(3.)*etav2_x7-etav2_z7)/2.*(2*is-1)*nn2_ijk(i,2)*(0,1)
endif

if (i>8) then	!101
phi_sigma_alpha2(0,0,is,is,i)=(+etav2_x1+sqrt(3.)*etav2_x2+sqrt(3.)*etav2_z1+3*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(0,1,is,is,i)=(-sqrt(3.)*etav2_x1+etav2_x2-3*etav2_z1+sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,is,i)=(-etav2_x7-sqrt(3.)*etav2_z7)/2.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(1,0,is,is,i)=(-sqrt(3.)*etav2_x1-3*etav2_x2+etav2_z1+sqrt(3.)*etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
phi_sigma_alpha2(1,1,is,is,i)=(3*etav2_x1-sqrt(3.)*etav2_x2-sqrt(3.)*etav2_z1+etav2_z2)/4.*(2*is-1)*nn2_ijk(i,3)
if (l_orbital_f>2) phi_sigma_alpha2(1,2,is,is,i)=(sqrt(3.)*etav2_x7-etav2_z7)/2.*(2*is-1)*nn2_ijk(i,3)

phi_sigma_alpha2(0,0,is,1-is,i)=(+etav2_x1-sqrt(3.)*etav2_x2-sqrt(3.)*etav2_z1+3*etav2_z2)/4.*nn2_ijk(i,1)
phi_sigma_alpha2(0,1,is,1-is,i)=(-sqrt(3.)*etav2_x1-etav2_x2+3*etav2_z1+sqrt(3.)*etav2_z2)/4.*nn2_ijk(i,1)
if (l_orbital_f>2) phi_sigma_alpha2(0,2,is,1-is,i)=(-etav2_x7+sqrt(3.)*etav2_z7)/2.*nn2_ijk(i,1)
phi_sigma_alpha2(1,0,is,1-is,i)=(-sqrt(3.)*etav2_x1+3*etav2_x2-etav2_z1+sqrt(3.)*etav2_z2)/4.*nn2_ijk(i,1)
phi_sigma_alpha2(1,1,is,1-is,i)=(+3*etav2_x1+sqrt(3.)*etav2_x2+sqrt(3.)*etav2_z1+etav2_z2)/4.*nn2_ijk(i,1)
if (l_orbital_f>2) phi_sigma_alpha2(1,2,is,1-is,i)=(+sqrt(3.)*etav2_x7+etav2_z7)/2.*nn2_ijk(i,1)

endif

enddo

do o1=0,l_orbital-1
do o2=0,l_orbital_f-1
do is=0,1
do is2=0,1
phi_alpha_sigma2(o2,o1,is2,is,i)=-conjg(phi_sigma_alpha2(o1,o2,is,is2,i))
enddo
enddo
enddo
enddo

enddo
endif


if(phi_type=="sm2") then	
  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
!xy !4*s_x sigma_x (c_y+c_z) and cyclic --> no along G-R and X-M
if(i<5)    phi_sigma_alpha2(0,0,is,1-is,i)=(nn2_ijk(i,1)-(0,1)*(2*is-1)*nn2_ijk(i,2))*eta_v2
!yz
if(i>4 .and. i<9)    phi_sigma_alpha2(0,0,is,1-is,i)=(-(0,1)*(2*is-1)*nn2_ijk(i,2))*eta_v2
if(i>4 .and. i<9)    phi_sigma_alpha2(0,0,is,is,i)=((2*is-1)*nn2_ijk(i,3))*eta_v2
!zx
if(i>8)    phi_sigma_alpha2(0,0,is,1-is,i)=nn2_ijk(i,1)*eta_v2
if(i>8)    phi_sigma_alpha2(0,0,is,is,i)=nn2_ijk(i,3)*(2*is-1)*eta_v2

!xy !s_x sigma_x c_y and cyclic 
!if(i<5)    phi_sigma_alpha2(0,0,is,1-is,i)=(nn2_ijk(i,1))*eta_v2
!yz
!if(i>4 .and. i<9)    phi_sigma_alpha2(0,0,is,1-is,i)=(-(0,1)*(2*is-1)*nn2_ijk(i,2))*eta_v2
!zx
!if(i>8)    phi_sigma_alpha2(0,0,is,is,i)=nn2_ijk(i,3)*(2*is-1)*eta_v2

!xy !s_x sigma_x (c_y-c_z) and cyclic --> no along G-R
!if(i<5)    phi_sigma_alpha2(0,0,is,1-is,i)=(nn2_ijk(i,1)+(0,1)*(2*is-1)*nn2_ijk(i,2))*eta_v2
!yz
!if(i>4 .and. i<9)    phi_sigma_alpha2(0,0,is,1-is,i)=-((0,1)*(2*is-1)*nn2_ijk(i,2))*eta_v2
!if(i>4 .and. i<9)    phi_sigma_alpha2(0,0,is,is,i)=-((2*is-1)*nn2_ijk(i,3))*eta_v2
!zx
!if(i>8)    phi_sigma_alpha2(0,0,is,1-is,i)=-nn2_ijk(i,1)*eta_v2
!if(i>8)    phi_sigma_alpha2(0,0,is,is,i)=nn2_ijk(i,3)*(2*is-1)*eta_v2


enddo

 do is=0,1	!sigma	!alpha
  do is2=0,1	!sigma	!alpha
phi_alpha_sigma2(0,0,is2,is,i)=-conjg(phi_sigma_alpha2(0,0,is,is2,i))
enddo
enddo

enddo
endif


!!!!!!!!!!!3rd nn



phi_sigma_alpha3=0
phi_alpha_sigma3=0


  do i=1, n_nn3
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0, l_orbital-1		!d
      do o2=0, l_orbital-1		!f for now only gamma8
    
     if (lphioverlap) then
      do k=-3,3
       phi_sigma_alpha3(0,0,is,is2,i)=phi_sigma_alpha3(0,0,is,is2,i)+alpha_msigma(0,is2,is,k)*sph_arm(j1,k,nn3_theta_phi(i,1),nn3_theta_phi(i,2))*sqrt(pi/3)*eta_v3
       phi_alpha_sigma3(0,0,is,is2,i)=phi_alpha_sigma3(0,0,is,is2,i)+alpha_msigma(0,is,is2,k)*conjg(sph_arm(j1,k,pi-nn3_theta_phi(i,1),pi+nn3_theta_phi(i,2)))*sqrt(pi/3)*eta_v3
      enddo
endif
enddo
enddo
enddo
enddo
enddo


if (l_orbital==2 .and. l_orbital_f==3) then
  do i=1, n_nn3
   do is=0,1	!sigma	!alpha

phi_sigma_alpha3(0,2,is,is,i)=2*etav3_x7*(2*is-1)*nn3_ijk(i,3)
phi_sigma_alpha3(0,2,is,1-is,i)=-etav3_x7*nn3_ijk(i,1)
phi_sigma_alpha3(0,2,is,1-is,i)=phi_sigma_alpha3(0,2,is,1-is,i)+etav3_x7*(0,1)*(2*is-1)*nn3_ijk(i,2)

phi_sigma_alpha3(1,2,is,1-is,i)=-etav3_z7*nn3_ijk(i,1)
phi_sigma_alpha3(1,2,is,1-is,i)=phi_sigma_alpha3(1,2,is,1-is,i)-etav3_z7*(0,1)*(2*is-1)*nn3_ijk(i,2)


phi_sigma_alpha3(0,0,is,is,i)=-etav3_x1*(2*is-1)*nn3_ijk(i,3)
phi_sigma_alpha3(0,0,is,1-is,i)=-etav3_x1b*nn3_ijk(i,1)
phi_sigma_alpha3(0,0,is,1-is,i)=phi_sigma_alpha3(0,0,is,1-is,i)+(0,1)*etav3_x1b*nn3_ijk(i,2)*(2*is-1)
phi_sigma_alpha3(1,1,is,is,i)=-etav3_z2*(2*is-1)*nn3_ijk(i,3)
phi_sigma_alpha3(1,1,is,1-is,i)=-etav3_z2b*nn3_ijk(i,1)
phi_sigma_alpha3(1,1,is,1-is,i)=phi_sigma_alpha3(1,1,is,1-is,i)+(0,1)*nn3_ijk(i,2)*etav3_z2b*(2*is-1)
phi_sigma_alpha3(0,1,is,is,i)=etav3_x2*(0,1)*nn3_ijk(i,1)*nn3_ijk(i,2)*nn3_ijk(i,3)
phi_sigma_alpha3(0,1,is,1-is,i)=etav3_x2b*nn3_ijk(i,1)
phi_sigma_alpha3(0,1,is,1-is,i)=phi_sigma_alpha3(0,1,is,1-is,i)+(0,1)*etav3_x2b*nn3_ijk(i,2)*(2*is-1)
phi_sigma_alpha3(1,0,is,is,i)=-etav3_x2*(0,1)*nn3_ijk(i,1)*nn3_ijk(i,2)*nn3_ijk(i,3)
phi_sigma_alpha3(1,0,is,1-is,i)=etav3_x2b*nn3_ijk(i,1)
phi_sigma_alpha3(1,0,is,1-is,i)=phi_sigma_alpha3(1,0,is,1-is,i)+(0,1)*etav3_x2b*nn3_ijk(i,2)*(2*is-1)

enddo

do o1=0,1
do o2=0,2
do is=0,1
do is2=0,1
phi_alpha_sigma3(o2,o1,is2,is,i)=-conjg(phi_sigma_alpha3(o1,o2,is,is2,i))
enddo
enddo
enddo
enddo

enddo





endif

if (phi_type=="pd2" .or. phi_type=="pd4") then

phi_sigma_alpha=0
phi_sigma_alpha2=0
phi_sigma_alpha3=0
phi_alpha_sigma=0
phi_alpha_sigma2=0
phi_alpha_sigma3=0

  do i=1, n_nn
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0, l_orbital-1		!d
      do o2=0, l_orbital-1		!f for now only gamma8
    
     if (lphioverlap) then
      do k=-3,3
 !      phi_sigma_alpha(0,0,is,is2,i)=phi_sigma_alpha(0,0,is,is2,i)+alpha_msigma(0,is2,is,k)*sph_arm(j1,k,nn_theta_phi(i,1),nn_theta_phi(i,2))*eta_v*&
  !     									sph_arm_d(1,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2))
   !    phi_alpha_sigma(0,0,is,is2,i)=phi_alpha_sigma(0,0,is,is2,i)+alpha_msigma(0,is,is2,k)*conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*eta_v*&
    !   									sph_arm_d(1,nn_theta_phi(i,1),nn_theta_phi(i,2))


       phi_sigma_alpha(o1,o2,is,is2,i)=phi_sigma_alpha(o1,o2,is,is2,i)+alpha_msigma(o2,is2,is,k)*sph_arm(j1,k,nn_theta_phi(i,1),nn_theta_phi(i,2))*&
       				sph_arm_d(1,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2))*eta_v
       phi_alpha_sigma(o1,o2,is,is2,i)=phi_alpha_sigma(o1,o2,is,is2,i)+alpha_msigma(o1,is,is2,k)*&
       				conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*sph_arm_d(1,nn_theta_phi(i,1),nn_theta_phi(i,2))*eta_v

       									
      enddo
endif
enddo
enddo
enddo
enddo
enddo


  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0, l_orbital-1		!d
      do o2=0, l_orbital-1		!f for now only gamma8
    
     if (lphioverlap) then
      do k=-3,3
       phi_sigma_alpha2(0,0,is,is2,i)=phi_sigma_alpha2(0,0,is,is2,i)+alpha_msigma(0,is2,is,k)*sph_arm(j1,k,nn2_theta_phi(i,1),nn2_theta_phi(i,2))*eta_v2*&
       									sph_arm_d(1,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2))
       phi_alpha_sigma2(0,0,is,is2,i)=phi_alpha_sigma2(0,0,is,is2,i)+alpha_msigma(0,is,is2,k)*conjg(sph_arm(j1,k,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2)))*eta_v2*&
       									sph_arm_d(1,nn2_theta_phi(i,1),nn2_theta_phi(i,2))
      enddo
endif
enddo
enddo
enddo
enddo
enddo




endif

if(phi_type=="sm2") then	
  do i=1, n_nn3
   do is=0,1	!sigma	!alpha

phi_sigma_alpha3(0,0,is,1-is,i)=nn3_ijk(i,1)*nn3_ijk(i,2)*nn3_ijk(i,3)*(1+(0,1)*(2*is-1))*eta_v3	!sx sy sz (sigma_x + sigma_y + sigma_z)
phi_sigma_alpha3(0,0,is,is,i)=nn3_ijk(i,1)*nn3_ijk(i,2)*nn3_ijk(i,3)*(2*is-1)*eta_v3


enddo

 do is=0,1	!sigma	!alpha
  do is2=0,1	!sigma	!alpha
phi_alpha_sigma3(0,0,is2,is,i)=-conjg(phi_sigma_alpha3(0,0,is,is2,i))
enddo
enddo

enddo
endif



end subroutine set_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kin_en !(kin_energy)
!complex(KIND=idpc):: kin_energy(0:l_orbital-1,0:l_orbital-1,0:1,0:1,n_nn,0:1)
integer		::ikx,iky,ikz


kin_energy=0
kin_energy2=0
kin_energy3=0


if (l_orbital==2) then

  do i=1, n_nn
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
    
   do o1=0,l_orbital_f-1	!orb 1
    do o2=0,l_orbital_f-1	!orb 2

!eta=from dzero paper
!eta1=multiplies the other term
!I divide by 4 at the end!


if( l_k_dz) then	!from prl cubic top ins
!x
if (is==is2 .and. i<3 .and. o1==0 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=1*eta1_df(:)+3*eta_df(:)
if (is==is2 .and. i<3 .and. o1==1 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=3*eta1_df(:)+eta_df(:)
if (is==is2 .and. i<3 .and. o1==0 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=sqrt(3.)*(1*eta1_df(:)-eta_df(:))
if (is==is2 .and. i<3 .and. o1==1 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=sqrt(3.)*(1*eta1_df(:)-eta_df(:))
!y
if (is==is2 .and. i<5 .and. i>2 .and. o1==0 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=1*eta1_df(:)+3*eta_df(:)
if (is==is2 .and. i<5 .and. i>2 .and. o1==1 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=3*eta1_df(:)+eta_df(:)
if (is==is2 .and. i<5 .and. i>2 .and. o1==0 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=-sqrt(3.)*(1*eta1_df(:)-eta_df(:))
if (is==is2 .and. i<5 .and. i>2 .and. o1==1 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=-sqrt(3.)*(1*eta1_df(:)-eta_df(:))
!z
if (is==is2 .and. i>4 .and. o1==0 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=4*eta1_df(:)	!x^2-y^2 eta1: etax, etaf1
if (is==is2 .and. i>4 .and. o1==1 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=4*eta_df(:)	!z^2	 eta : etaz, etaf2

else	!from overlap of wfc

!d
if (is==is2)   kin_energy(o1,o2,is,is2,i,0)=kin_energy(o1,o2,is,is2,i,0)+sph_arm_d(o1,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2))*sph_arm_d(o2,nn_theta_phi(i,1),nn_theta_phi(i,2))*&
						eta_df(0)*64*pi/5	!d

!f
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy(o1,o2,is,is2,i,1)=kin_energy(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn_theta_phi(i,1),nn_theta_phi(i,2)))*&
   							     eta_df(1)*16*pi/3		!f, then I divide by 4


       enddo
      enddo
     enddo

endif	!l_k_dz

enddo	!is
enddo	!is2
enddo	!o1
enddo	!o2

if(phi_type=="sm7" .or. phi_type=="smb") then

   do is=0,1	!sigma	!alpha
is2=is    
   do o1=0,l_orbital_f-1	!orb 1
    do o2=0,l_orbital_f-1	!orb 2

!Gamma7
!with itself
if (i<3 .and. o1==2 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=4*eta_7
if (i<5 .and. i>2  .and. o1==2 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=4*eta_7
if (i>4 .and. o1==2 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=4*eta_7

if (.not. l7z) then
!with Gamma8(1)
!x
if (i<3 .and. o1==0 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<3 .and. o1==1 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*sqrt(3.)*eta_78
if (i<3 .and. o1==2 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<3 .and. o1==2 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=-2*sqrt(3.)*eta_78
!y
if (i<5 .and. i>2  .and. o1==0 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<5 .and. i>2  .and. o1==1 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=+2*sqrt(3.)*eta_78
if (i<5 .and. i>2  .and. o1==2 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<5 .and. i>2  .and. o1==2 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=+2*sqrt(3.)*eta_78
!z
if (i>4 .and. o1==0 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=+4*eta_78
if (i>4 .and. o1==2 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=+4*eta_78

else

!with Gamma8(2)
!x
if (i<3 .and. o1==1 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<3 .and. o1==0 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=+2*sqrt(3.)*eta_78
if (i<3 .and. o1==2 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<3 .and. o1==2 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=+2*sqrt(3.)*eta_78
!y
if (i<5 .and. i>2  .and. o1==1 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<5 .and. i>2  .and. o1==0 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=-2*sqrt(3.)*eta_78
if (i<5 .and. i>2  .and. o1==2 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=-2*eta_78
if (i<5 .and. i>2  .and. o1==2 .and. o2==0)   kin_energy(o1,o2,is,is2,i,:)=-2*sqrt(3.)*eta_78
!z
if (i>4 .and. o1==1 .and. o2==2)   kin_energy(o1,o2,is,is2,i,:)=+4*eta_78
if (i>4 .and. o1==2 .and. o2==1)   kin_energy(o1,o2,is,is2,i,:)=+4*eta_78
endif



enddo
enddo
enddo

endif

enddo	!i

kin_energy=kin_energy/4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
    
   do o1=0,l_orbital-1	!orb 1
    do o2=0,l_orbital-1	!orb 2


!d
!2nd nearest neighbor: from wfc overlap
if (.not. l_k_dz) then
if (is==is2)   kin_energy2(o1,o2,is,is2,i,0)=kin_energy2(o1,o2,is,is2,i,0)+eta2_df(0)*sph_arm_d(o1,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2))*sph_arm_d(o2,nn2_theta_phi(i,1),nn2_theta_phi(i,2))&
									*4*16*pi/5	!d
else

!110 1234
!2nd nearest neighbor: analytical
if(i<5) kin_energy2(1,1,is,is,i,0)=eta2_df(0)!/4.
if(i<5) kin_energy2(0,0,is,is,i,0)=eta2_dx!/4.

!011 5678
if(i>4 .and. i<9) kin_energy2(0,0,is,is,i,0)=(3*eta2_df(0)+eta2_dx)/4.		!x2-y2
if(i>4 .and. i<9) kin_energy2(1,1,is,is,i,0)=(eta2_df(0)+3*eta2_dx)/4.		!z^2
if(i>4 .and. i<9) kin_energy2(0,1,is,is,i,0)=(eta2_dx-eta2_df(0))*sqrt(3.)/4.
if(i>4 .and. i<9) kin_energy2(1,0,is,is,i,0)=(eta2_dx-eta2_df(0))*sqrt(3.)/4.

!101 9101112
if(i>8) kin_energy2(0,0,is,is,i,0)=(3*eta2_df(0)+eta2_dx)/4.
if(i>8) kin_energy2(1,1,is,is,i,0)=(eta2_df(0)+3*eta2_dx)/4.
if(i>8) kin_energy2(0,1,is,is,i,0)=(eta2_df(0)-eta2_dx)*sqrt(3.)/4.
if(i>8) kin_energy2(1,0,is,is,i,0)=(eta2_df(0)-eta2_dx)*sqrt(3.)/4.
endif !l_k_dz


enddo	!o1
enddo	!o2

!f
if (.not. l_k_dz) then
   do o1=0,l_orbital_f-1	!orb 1
    do o2=0,l_orbital_f-1	!orb 2
if(o1==2 .and. o2==2) then
eta2_ok=eta2_f7
else if (o1<2 .and. o2==2) then
eta2_ok=eta2_f78
else if (o1==2 .and. o2<2) then
eta2_ok=eta2_f78
else
eta2_ok=eta2_f
endif

      do k=-3,3
       do k2=-3,3
        do is3=0,1
        
   kin_energy2(o1,o2,is,is2,i,1)=kin_energy2(o1,o2,is,is2,i,1)+eta2_ok*alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn2_theta_phi(i,1),nn2_theta_phi(i,2)))&		!f
   							     *4.0d0*pi/3.0d0

       enddo
      enddo
     enddo

    enddo	!o1
   enddo	!o2



else

!110 1234
if (i<5) then
kin_energy2(0,0,is,is,i,1)=eta2_f81
kin_energy2(1,1,is,is,i,1)=eta2_f82
if (l_orbital_f >2) kin_energy2(2,2,is,is,i,1)=eta2_f7
kin_energy2(0,1,is,is,i,1)=eta2_f812*(0,1)*(-1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,2)
kin_energy2(1,0,is,is,i,1)=eta2_f812*(0,1)*(+1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,2)
if (l_orbital_f >2) kin_energy2(0,2,is,is,i,1)=eta2_f781
if (l_orbital_f >2) kin_energy2(2,0,is,is,i,1)=eta2_f781
if (l_orbital_f >2) kin_energy2(1,2,is,is,i,1)=eta2_f782*(0,1)*(+1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,2)
if (l_orbital_f >2) kin_energy2(2,1,is,is,i,1)=eta2_f782*(0,1)*(-1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,2)
endif

!011 5678
if(i>4 .and. i<9) then
kin_energy2(0,0,is,is,i,1)=(eta2_f81+3*eta2_f82)/4.
kin_energy2(1,1,is,is,i,1)=(3*eta2_f81+eta2_f82)/4.
if (l_orbital_f >2) kin_energy2(2,2,is,is,i,1)=eta2_f7
kin_energy2(0,1,is,is,i,1)=sqrt(3.)/4.*(eta2_f81-eta2_f82)
kin_energy2(1,0,is,is,i,1)=sqrt(3.)/4.*(eta2_f81-eta2_f82)
kin_energy2(0,1,is,1-is,i,1)=eta2_f812*(0,1)*(-1)*nn2_ijk(i,2)*nn2_ijk(i,3)
kin_energy2(1,0,is,1-is,i,1)=eta2_f812*(0,1)*(+1)*nn2_ijk(i,2)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(0,2,is,is,i,1)=-eta2_f781/2.
if (l_orbital_f >2) kin_energy2(2,0,is,is,i,1)=-eta2_f781/2.
if (l_orbital_f >2) kin_energy2(0,2,is,1-is,i,1)=eta2_f782*(0,1)*(+1)*sqrt(3.)/2.*nn2_ijk(i,2)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(2,0,is,1-is,i,1)=eta2_f782*(0,1)*(-1)*sqrt(3.)/2.*nn2_ijk(i,2)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(1,2,is,is,i,1)=-sqrt(3.)/2.*eta2_f781
if (l_orbital_f >2) kin_energy2(2,1,is,is,i,1)=-sqrt(3.)/2.*eta2_f781
if (l_orbital_f >2) kin_energy2(1,2,is,1-is,i,1)=eta2_f782*(0,1)*(-1)/2.*nn2_ijk(i,2)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(2,1,is,1-is,i,1)=eta2_f782*(0,1)*(+1)/2.*nn2_ijk(i,2)*nn2_ijk(i,3)
endif
!101 9101112
if(i>8) then
kin_energy2(0,0,is,is,i,1)=(eta2_f81+3*eta2_f82)/4.
kin_energy2(1,1,is,is,i,1)=(3*eta2_f81+eta2_f82)/4.
if (l_orbital_f >2) kin_energy2(2,2,is,is,i,1)=eta2_f7
kin_energy2(0,1,is,is,i,1)=-sqrt(3.)/4.*(eta2_f81-eta2_f82)
kin_energy2(1,0,is,is,i,1)=-sqrt(3.)/4.*(eta2_f81-eta2_f82)
kin_energy2(0,1,is,1-is,i,1)=eta2_f812*(-1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,3)
kin_energy2(1,0,is,1-is,i,1)=eta2_f812*(+1)*(2*is-1)*nn2_ijk(i,1)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(0,2,is,is,i,1)=-eta2_f781/2.
if (l_orbital_f >2) kin_energy2(2,0,is,is,i,1)=-eta2_f781/2.
if (l_orbital_f >2) kin_energy2(0,2,is,1-is,i,1)=eta2_f782*(-1)*(2*is-1)*sqrt(3.)/2.*nn2_ijk(i,1)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(2,0,is,1-is,i,1)=eta2_f782*(+1)*(2*is-1)*sqrt(3.)/2.*nn2_ijk(i,1)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(1,2,is,is,i,1)=sqrt(3.)/2.*eta2_f781
if (l_orbital_f >2) kin_energy2(2,1,is,is,i,1)=sqrt(3.)/2.*eta2_f781
if (l_orbital_f >2) kin_energy2(1,2,is,1-is,i,1)=eta2_f782*(-1)*(2*is-1)/2.*nn2_ijk(i,1)*nn2_ijk(i,3)
if (l_orbital_f >2) kin_energy2(2,1,is,1-is,i,1)=eta2_f782*(+1)*(2*is-1)/2.*nn2_ijk(i,1)*nn2_ijk(i,3)



endif

endif !l_k_dz


enddo !is2
enddo !is
enddo !i


!print *, kin_energy2(2,2,0,0,1,1), kin_energy2(2,2,0,0,5,1)
!print *, kin_energy2(2,2,0,1,1,1), kin_energy2(2,2,0,1,5,1)
!print *, kin_energy2(2,2,1,0,1,1), kin_energy2(2,2,1,0,5,1)
!print *, kin_energy2(2,2,1,1,1,1), kin_energy2(2,2,1,1,5,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!3rd

kin_energy3=0


  do i=1, n_nn3

!f
if (.not. l_k_dz) then
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
     do o1=0,l_orbital_f-1	!orb 1
      do o2=0,l_orbital_f-1	!orb 2
    
if(o1==2 .and. o2==2) then
eta3_ok=eta3_f7*(9.d0/5.0d0)
else if (o1<2 .and. o2==2) then
eta3_ok=eta3_78*(9.0/sqrt(10.0d0))		!eta3_78
else if (o1==2 .and. o2<2) then
eta3_ok=eta3_78*(9.0/sqrt(10.0d0))		!eta3_78
else
eta3_ok=eta3_f8*(9.0d0/2.0d0)	!to have 1
endif

      do k=-3,3
       do k2=-3,3
        do is3=0,1
        
   kin_energy3(o1,o2,is,is2,i,1)=kin_energy3(o1,o2,is,is2,i,1)-eta3_ok*alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn3_theta_phi(i,1),pi+nn3_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn3_theta_phi(i,1),nn3_theta_phi(i,2)))&		!f
   							     *4.0d0*pi/3.0d0


       enddo
      enddo
     enddo
     
    enddo
   enddo
  enddo
 enddo

else

   do is=0,1	!sigma	!alpha
!    do is2=0,1	!alpha	!sigma
!     do o1=0,l_orbital_f-1	!orb 1
!      do o2=0,l_orbital_f-1	!orb 2

kin_energy3(0,0,is,is,i,1)=eta3_f8
kin_energy3(1,1,is,is,i,1)=eta3_f8
if (l_orbital_f >2) kin_energy3(2,2,is,is,i,1)=eta3_f7

!I ignore nondiagonal in Gamma8 sector (it is zero from dft) (?)
!kin_energy3(0,1,is,is,i,1)=+eta3_f8*(0,1)/sqrt(3.)*(2*is-1)*nn3_ijk(i,1)*nn3_ijk(i,2)
!kin_energy3(1,0,is,is,i,1)=-eta3_f8*(0,1)/sqrt(3.)*(2*is-1)*nn3_ijk(i,1)*nn3_ijk(i,2)
!kin_energy3(0,1,is,1-is,i,1)=-eta3_f8/sqrt(3.)*(2*is-1)*nn3_ijk(i,1)*nn3_ijk(i,3)
!kin_energy3(1,0,is,1-is,i,1)=+eta3_f8/sqrt(3.)*(2*is-1)*nn3_ijk(i,1)*nn3_ijk(i,3)
!kin_energy3(0,1,is,1-is,i,1)=kin_energy3(0,1,is,1-is,i,1)+eta3_f8*(0,1)/sqrt(3.)*nn3_ijk(i,2)*nn3_ijk(i,3)
!kin_energy3(1,0,is,1-is,i,1)=kin_energy3(1,0,is,1-is,i,1)-eta3_f8*(0,1)/sqrt(3.)*nn3_ijk(i,2)*nn3_ijk(i,3)


!Nondiagonal Gamma8-Gamma7
if (l_orbital_f >2) kin_energy3(0,2,is,1-is,i,1)=eta3_78*(-nn3_ijk(i,1)*nn3_ijk(i,3)*(2*is-1)+(0,1)*nn3_ijk(i,2)*nn3_ijk(i,3))
if (l_orbital_f >2) kin_energy3(2,0,is,1-is,i,1)=eta3_78*(+nn3_ijk(i,1)*nn3_ijk(i,3)*(2*is-1)-(0,1)*nn3_ijk(i,2)*nn3_ijk(i,3))

if (l_orbital_f >2) kin_energy3(1,2,is,is,i,1)=2*eta3_78/sqrt(3.)*(0,1)*nn3_ijk(i,1)*nn3_ijk(i,2)*(2*is-1)
if (l_orbital_f >2) kin_energy3(1,2,is,1-is,i,1)=eta3_78/sqrt(3.)*(-nn3_ijk(i,1)*nn3_ijk(i,3)*(2*is-1)-(0,1)*nn3_ijk(i,2)*nn3_ijk(i,3))
if (l_orbital_f >2) kin_energy3(2,1,is,is,i,1)=-2*eta3_78/sqrt(3.)*(0,1)*nn3_ijk(i,1)*nn3_ijk(i,2)*(2*is-1)
if (l_orbital_f >2) kin_energy3(2,1,1-is,is,i,1)=eta3_78/sqrt(3.)*(-nn3_ijk(i,1)*nn3_ijk(i,3)*(2*is-1)+(0,1)*nn3_ijk(i,2)*nn3_ijk(i,3))

!    enddo
!   enddo
!  enddo
 enddo



endif !lkdz

!d

  do is=0,1	!sigma	!alpha
   do o1=0,l_orbital-1	!orb 1
   kin_energy3(o1,o1,is,is,i,0)=eta3_d
  enddo
 enddo

enddo	!i

!print *, kin_energy3(2,2,0,0,1,1)



!only for Gamma7
!2nd
!if (l_orbital_f==3) then
!do i=1, n_nn2
!kin_energy2(2,2,0,0,i,1)=eta2_f7
!kin_energy2(2,2,1,1,i,1)=eta2_f7
!kin_energy2(2,2,1,0,i,1)=0
!kin_energy2(2,2,0,1,i,1)=0
!enddo
!3rd
!do i=1, n_nn3
!kin_energy3(2,2,0,0,i,1)=eta3_f7
!kin_energy3(2,2,1,1,i,1)=eta3_f7
!kin_energy3(2,2,1,0,i,1)=0
!kin_energy3(2,2,0,1,i,1)=0
!enddo
!endif

!!!!!!!!!!!!!!!!!!!1 orb model

!l_orbital==1
else
!kin_energy2=0
!kin_energy3=0
do is=0,1
!kin_energy(0,0,is,is,:,:)=eta_7
!kin_energy2(0,0,is,is,:,:)=eta2_f7
!kin_energy3(0,0,is,is,:,:)=eta3_f7
do i=1, n_nn
kin_energy(0,0,is,is,i,:)=eta1_d
!if (i>4) kin_energy(0,0,is,is,i,1)=4
!if (i<5) kin_energy(0,0,is,is,i,:)=1
!if (i>4) kin_energy(0,0,is,is,i,0)=1
!if (i>4) kin_energy(0,0,is,is,i,1)=4
enddo
do i=1, n_nn2
kin_energy2(0,0,is,is,i,:)=eta2_d
enddo
do i=1, n_nn3
kin_energy3(0,0,is,is,i,:)=eta3_d
enddo
enddo

!2nd
if (phi_type=="g7_" .or. phi_type=="sm2" .or. phi_type=="trv" .or. phi_type=="ons") then
do i=1, n_nn
kin_energy(0,0,0,0,i,1)=eta_7
kin_energy(0,0,1,1,i,1)=eta_7
enddo
do i=1, n_nn2
kin_energy2(0,0,0,0,i,1)=eta2_f7
kin_energy2(0,0,1,1,i,1)=eta2_f7
enddo
!3rd
do i=1, n_nn3
kin_energy3(0,0,0,0,i,1)=eta3_f7
kin_energy3(0,0,1,1,i,1)=eta3_f7
enddo
endif


!!!kin en by overlap
o1=0
o2=0
if (lkinoverlap .and. .not. phi_type=="g7_") then

kin_energy(:,:,:,:,:,1)=0
  do i=1, n_nn
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy(0,0,is,is2,i,1)=kin_energy(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn_theta_phi(i,1),nn_theta_phi(i,2)))*&
   							     eta_7*4*pi/3		


       enddo
      enddo
     enddo
enddo
enddo
enddo

kin_energy2(:,:,:,:,:,1)=0
  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy2(0,0,is,is2,i,1)=kin_energy2(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn2_theta_phi(i,1),nn2_theta_phi(i,2)))*&
   							     eta2_f7*4*pi/3		


       enddo
      enddo
     enddo
      enddo
     enddo
enddo

kin_energy3(:,:,:,:,:,1)=0
  do i=1, n_nn3
   do is=0,1	!sigma	!alpha
    do is2=0,1	!alpha	!sigma
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy3(o1,o2,is,is2,i,1)=kin_energy3(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn3_theta_phi(i,1),pi+nn3_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn3_theta_phi(i,1),nn3_theta_phi(i,2)))*&
   							     eta3_f7*4*pi/3		


       enddo
      enddo
      enddo
     enddo
     enddo
enddo

endif


endif

if (phi_type=="pd2" .or. phi_type=="pd4") then
kin_energy=0
kin_energy2=0
kin_energy3=0
o1=0
o2=0

!d1
  do i=1, n_nn
   do is=0,1	!sigma	!alpha
kin_energy(o1,o2,is,is,i,0)=kin_energy(o1,o2,is,is,i,0)+eta1_d*sph_arm_d(1,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2))*sph_arm_d(1,nn_theta_phi(i,1),nn_theta_phi(i,2))!*&
							!d
enddo
enddo

!d2
  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
kin_energy2(o1,o2,is,is,i,0)=kin_energy2(o1,o2,is,is,i,0)+eta2_d*sph_arm_d(1,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2))*sph_arm_d(1,nn2_theta_phi(i,1),nn2_theta_phi(i,2))!&
										!d
enddo
enddo

!f1
  do i=1, n_nn
   do is=0,1	!sigma	!alpha
   do is2=0,1	!sigma	!alpha
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy(o1,o2,is,is2,i,1)=kin_energy(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn_theta_phi(i,1),pi+nn_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn_theta_phi(i,1),nn_theta_phi(i,2)))*&
   							     eta_7*16*pi/3		!f, then I divide by 4


       enddo
      enddo
     enddo
       enddo
      enddo
     enddo

!f2
  do i=1, n_nn2
   do is=0,1	!sigma	!alpha
   do is2=0,1	!sigma	!alpha
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy2(o1,o2,is,is2,i,1)=kin_energy2(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn2_theta_phi(i,1),pi+nn2_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn2_theta_phi(i,1),nn2_theta_phi(i,2)))*&
   							     eta2_f7*16*pi/3		!f, then I divide by 4


       enddo
      enddo
     enddo
       enddo
      enddo
     enddo


!f3
  do i=1, n_nn3
   do is=0,1	!sigma	!alpha
   do is2=0,1	!sigma	!alpha
      do k=-3,3
       do k2=-3,3
        do is3=0,1
   kin_energy3(o1,o2,is,is2,i,1)=kin_energy3(o1,o2,is,is2,i,1)+alpha_msigma(o1,is,is3,k)*  conjg(sph_arm(j1,k,pi-nn3_theta_phi(i,1),pi+nn3_theta_phi(i,2)))*&
   							     alpha_msigma(o2,is2,is3,k2)*(sph_arm(j1,k2,nn3_theta_phi(i,1),nn3_theta_phi(i,2)))*&
   							     eta3_f7*16*pi/3		!f, then I divide by 4


       enddo
      enddo
     enddo
       enddo
      enddo
     enddo

endif

if (phi_type=="q2d") then
  do is=0,1
  do i=1, n_nn
if (i<5)   kin_energy(0,0,is,is,i,1)=eta_7	!100 010
if (i>4)   kin_energy(0,0,is,is,i,1)=0 !0.1*eta_7	!001
enddo
  do i=1, n_nn2
if (i<5)   kin_energy2(0,0,is,is,i,1)=eta2_f7	!110
if (i>4)   kin_energy2(0,0,is,is,i,1)=0 !0.1*eta2_f7!101 011
enddo

enddo
endif

end subroutine set_kin_en


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_phi_k(phi_k_sa,phi_k_as, ke_k, delta_k,nk,k_vecs_in,d23)
complex(KIND=idpc),intent(out)	:: phi_k_sa(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1),phi_k_as(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1),ke_k(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,0:1)
real(KIND=idp),intent(out)	:: delta_k(0:nk-1)
integer				:: i,is,is2,ik, nn_max, nn2_max, nn3_max
integer,intent(in)		:: nk,d23
real(kind=idp),intent(in)	:: k_vecs_in(0:nk-1,4)



if (d23==3) then
 nn_max=n_nn
 nn2_max=n_nn2
 nn3_max=n_nn3
else if (d23==2) then	!only when z=0
 if (lkz) nn_max=n_nn-2	!100 -100 010 0-10
 if (lkz) nn2_max=n_nn2-8!110 1-10 -110 -1-10
 if (lkz) nn3_max=0!noone in plane
 if (.not. lkz) nn_max=n_nn
 if (.not. lkz) nn2_max=n_nn2!?
 if (.not. lkz) nn3_max=n_nn3!?
else
 print *, "error set_phi_k"
 stop
endif

!print *, "set_phi_k", nn_max, d23

phi_k_sa=0
phi_k_as=0
ke_k=0
delta_k=0

do ik=0, nk-1

if (Vtki) then


!D=3
if (d23==3) then

!1st nn
 do is=0,1	!sigma !alpha
  do is2=0,1	!alpha	!sigma
  
   do i=1, nn_max

    do o1=0,l_orbital_f-1
     do o2=0,l_orbital_f-1

     ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)+nn_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     kin_energy(o1,o2,is,is2,i,:)
     enddo
    enddo

    do o1=0,l_orbital-1
     do o2=0,l_orbital_f-1

      phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)+nn_ijk(i,3)*k_vecs_in(ik,3)))*&
		    		     phi_sigma_alpha(o1,o2,is,is2,i)
      phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)+nn_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     phi_alpha_sigma(o2,o1,is,is2,i)
     enddo
    enddo
 
   enddo !i

!2nd nn
    do i=1, nn2_max
  
     do o1=0,l_orbital_f-1
      do o2=0,l_orbital_f-1

       ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)+nn2_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     kin_energy2(o1,o2,is,is2,i,:)

  
  
      enddo
     enddo


     do o1=0,l_orbital-1
      do o2=0,l_orbital_f-1

  
        phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)+nn2_ijk(i,3)*k_vecs_in(ik,3)))*&
		    		     phi_sigma_alpha2(o1,o2,is,is2,i)
        phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)+nn2_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     phi_alpha_sigma2(o2,o1,is,is2,i)
  
      enddo
     enddo
  
    enddo !i


!3rd nn
    do i=1, nn3_max

     do o1=0,l_orbital_f-1
      do o2=0,l_orbital_f-1

       ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)+nn3_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     kin_energy3(o1,o2,is,is2,i,:)

      enddo
     enddo

     do o1=0,l_orbital-1
      do o2=0,l_orbital_f-1

  
        phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)+nn3_ijk(i,3)*k_vecs_in(ik,3)))*&
		    		     phi_sigma_alpha3(o1,o2,is,is2,i)
        phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)+nn3_ijk(i,3)*k_vecs_in(ik,3)))*&
    				     phi_alpha_sigma3(o2,o1,is,is2,i)
  
      enddo
     enddo

    enddo !i


  enddo !is
 enddo	!is2

endif !D=3





!!!!!!!!!!D=2



if (d23==2) then
do is=0,1	!sigma !alpha
 do is2=0,1	!alpha	!sigma

!1st nn
  do i=1, nn_max

   do o1=0,l_orbital-1
    do o2=0, l_orbital_f-1

    phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)))*&
		    		     phi_sigma_alpha(o1,o2,is,is2,i)
    phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_alpha_sigma(o2,o1,is,is2,i)
    enddo
   enddo

   do o1=0,l_orbital_f-1
    do o2=0, l_orbital_f-1
     ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_in(ik,1)+nn_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     kin_energy(o1,o2,is,is2,i,:)
    enddo
   enddo

  enddo	!i

!2nd nn
 do i=1, nn2_max
 
  do o1=0,l_orbital_f-1
   do o2=0, l_orbital_f-1

     ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     kin_energy2(o1,o2,is,is2,i,:)

    enddo
   enddo

     do o1=0,l_orbital-1
      do o2=0,l_orbital_f-1

  
        phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
		    		     phi_sigma_alpha2(o1,o2,is,is2,i)
        phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_alpha_sigma2(o2,o1,is,is2,i)
  
      enddo
     enddo

    enddo !i

!3rd nn

 do i=1, nn3_max
 
  do o1=0,l_orbital_f-1
   do o2=0, l_orbital_f-1

     ke_k(o1,o2,is,is2,ik,:)=ke_k(o1,o2,is,is2,ik,:)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     kin_energy3(o1,o2,is,is2,i,:)

    enddo
   enddo

     do o1=0,l_orbital-1
      do o2=0,l_orbital_f-1

  
        phi_k_sa(o1,o2,is,is2,ik)=phi_k_sa(o1,o2,is,is2,ik)+&
    			    	     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
		    		     phi_sigma_alpha3(o1,o2,is,is2,i)
        phi_k_as(o2,o1,is,is2,ik)=phi_k_as(o2,o1,is,is2,ik)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_alpha_sigma3(o2,o1,is,is2,i)
  
      enddo
     enddo

    enddo !i




enddo !is2
enddo	!is
endif




else
phi_k_sa(0,0,0,0,ik)=1
phi_k_sa(0,0,1,1,ik)=1
phi_k_as(0,0,0,0,ik)=1
phi_k_as(0,0,1,1,ik)=1
delta_k(ik)=1
endif


enddo


end subroutine set_phi_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kin_energy2_kr(kin_energy2_kr,phi_kr_sa2,phi_kr_as2,  nk, k_vecs_in)
complex(kind=idpc) :: kin_energy2_kr(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn2,0:1), &
			phi_kr_sa2(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn2), phi_kr_as2(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1,n_nn2)
integer,intent(in)		:: nk
real(kind=idp),intent(in)	:: k_vecs_in(0:nk-1,4)


kin_energy2_kr=0
phi_kr_sa2=0
phi_kr_as2=0

do ik=0, nk-1
 do is=0,1	!sigma !alpha
  do is2=0,1	!alpha	!sigma
  do o1=0,l_orbital_f-1
   do o2=0,l_orbital_f-1
    do i=1, n_nn2	!nn

if (nn2_ijk(i,3)==1 .or. nn2_ijk(i,3)==-1)    kin_energy2_kr(o1,o2,is,is2,ik,i,:)=kin_energy2_kr(o1,o2,is,is2,ik,i,:)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     kin_energy2(o1,o2,is,is2,i,:)


if ((nn2_ijk(i,3)==1 .or. nn2_ijk(i,3)==-1) .and. o1<l_orbital)   &
				     phi_kr_sa2(o1,o2,is,is2,ik,i)=phi_kr_sa2(o1,o2,is,is2,ik,i)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_sigma_alpha2(o1,o2,is,is2,i)

if ((nn2_ijk(i,3)==1 .or. nn2_ijk(i,3)==-1)  .and. o2<l_orbital) &
				     phi_kr_as2(o1,o2,is,is2,ik,i)=phi_kr_as2(o1,o2,is,is2,ik,i)+&
    				     exp(2*(0,1)*pi*(nn2_ijk(i,1)*k_vecs_in(ik,1)+nn2_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_alpha_sigma2(o1,o2,is,is2,i)


enddo
enddo
enddo
enddo
enddo
enddo




end subroutine set_kin_energy2_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kin_energy3_kr(kin_energy3_kr,phi_kr_sa3,phi_kr_as3, nk, k_vecs_in)
complex(kind=idpc) :: kin_energy3_kr(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn3,0:1), &
			phi_kr_sa3(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn3), phi_kr_as3(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1,n_nn3)
integer,intent(in)		:: nk
real(kind=idp),intent(in)	:: k_vecs_in(0:nk-1,4)


kin_energy3_kr=0
phi_kr_sa3=0
phi_kr_as3=0

do ik=0, nk-1
 do is=0,1	!sigma !alpha
  do is2=0,1	!alpha	!sigma
  do o1=0,l_orbital_f-1
   do o2=0,l_orbital_f-1
    do i=1, n_nn3	!nn

 kin_energy3_kr(o1,o2,is,is2,ik,i,:)=kin_energy3_kr(o1,o2,is,is2,ik,i,:)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     kin_energy3(o1,o2,is,is2,i,:)
    				     
if ((nn3_ijk(i,3)==1 .or. nn3_ijk(i,3)==-1) .and. o1<l_orbital)   &
				     phi_kr_sa3(o1,o2,is,is2,ik,i)=phi_kr_sa3(o1,o2,is,is2,ik,i)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_sigma_alpha3(o1,o2,is,is2,i)

if ((nn3_ijk(i,3)==1 .or. nn3_ijk(i,3)==-1)  .and. o2<l_orbital) &
				     phi_kr_as3(o1,o2,is,is2,ik,i)=phi_kr_as3(o1,o2,is,is2,ik,i)+&
    				     exp(2*(0,1)*pi*(nn3_ijk(i,1)*k_vecs_in(ik,1)+nn3_ijk(i,2)*k_vecs_in(ik,2)))*&
    				     phi_alpha_sigma3(o1,o2,is,is2,i)

enddo
enddo
enddo
enddo
enddo
enddo




end subroutine set_kin_energy3_kr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,ldiag,nk,phi_sa, phi_as, ke, k_vecs)
complex(KIND=idpc),intent(out) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1)
complex(KIND=idpc),intent(in) 	:: phi_sa(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1),phi_as(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1),ke(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,0:1)
real(KIND=idp),intent(in) 	:: b_k, lambda_k, mu_k, k_vecs(0:nk-1,4)
integer, intent(in)		::ik,nk
logical,intent(in)		::ldiag

!order of states: c up, c down, f +, f -

!to redo with kinetic energy

Ham_k=0

!if (l_orbital==1) then
!if (ldiag) then
!Ham_k(0,0)=(ec_+tc_*factor_k(ik))-mu_k		!c
!Ham_k(1,1)=Ham_k(0,0)
!Ham_k(2,2)=(ef_+tf_*factor_k(ik)*b_k**2)-lambda_k-mu_k		!f
!if (.not. lmucf) Ham_k(2,2)=(ef_+tf_*factor_k(ik)*b_k**2)-lambda_k		!f
!Ham_k(3,3)=Ham_k(2,2)
!endif


!Ham_k(0,2)=Vcf_*phi_sa(0,0,0,0,ik)*b_k
!Ham_k(0,3)=Vcf_*phi_sa(0,0,0,1,ik)*b_k
!Ham_k(1,2)=Vcf_*phi_sa(0,0,1,0,ik)*b_k
!Ham_k(1,3)=Vcf_*phi_sa(0,0,1,1,ik)*b_k

!Ham_k(2,0)=Vcf_*phi_as(0,0,0,0,ik)*b_k
!Ham_k(2,1)=Vcf_*phi_as(0,0,0,1,ik)*b_k
!Ham_k(3,0)=Vcf_*phi_as(0,0,1,0,ik)*b_k
!Ham_k(3,1)=Vcf_*phi_as(0,0,1,1,ik)*b_k


!order of states: d1 down, d1 up, d2 down, d2 up,f1 -, f1 +, f2 -, f2 +, f3-, f3+

!elseif (l_orbital==2) then
if (.not. lreadham) then

do icf=0,1
do o1=0,l_orb(icf)-1
do is=0,1
ind1=index_cso(icf,is,o1)

do icf2=0,1
do o2=0,l_orb(icf2)-1
do is2=0,1

ind2=index_cso(icf2,is2,o2)

!kinetic energy
!if (o1==o2) Ham_k(ind1,ind2)=ke(o1,o2,is,is2,ik,o1)

if (ind1==ind2 .and. icf==0 .and. ldiag) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+ec_site_k(o1)-mu_k
if (ind1==ind2 .and. icf==1 .and. ldiag) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+ef_site_k(o1)-mu_k-lambda_k
if (icf==0 .and. icf2==0) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+tc_*ke(o1,o2,is,is2,ik,0)
if (icf==1 .and. icf2==1) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+b_k**2*tf_*ke(o1,o2,is,is2,ik,1)
if (icf==0 .and. icf2==1 .and. is==is2 .and. phi_type=="ons") Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+Vcf_*b_k
if (icf==1 .and. icf2==0 .and. is==is2 .and. phi_type=="ons") Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+Vcf_*b_k 


!sigma_alpha
if (icf==0 .and. icf2==1) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+phi_sa(o1,o2,is,is2,ik)*Vcf_*b_k

!alpha_sigma
if (icf==1 .and. icf2==0) Ham_k(ind1,ind2)=Ham_k(ind1,ind2)+phi_as(o1,o2,is,is2,ik)*Vcf_*b_k

enddo
enddo
enddo
enddo
enddo
enddo


else !lreadham

do lx=-n_shells,n_shells
do ly=-n_shells,n_shells
do lz=-n_shells,n_shells

if (sqrt(real(lx**2+ly**2+lz**2)) <=dist_max .and. sqrt(real(lx**2+ly**2+lz**2))>0.001) then

!dd
do i=0, 3
do j=0, 3
 Ham_k(i,j)=Ham_k(i,j)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly+k_vecs(ik,3)*lz)*(0,1))
enddo
enddo

!df
do i=0, 3
do j=4, dim_hk-1

 Ham_k(i,j)=Ham_k(i,j)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly+k_vecs(ik,3)*lz)*(0,1))*b_k
 Ham_k(j,i)=Ham_k(j,i)+Ham_input(lx,ly,lz,j,i)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly+k_vecs(ik,3)*lz)*(0,1))*b_k

enddo
enddo

!ff
do i=4, dim_hk-1
do j=4, dim_hk-1
 Ham_k(i,j)=Ham_k(i,j)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly+k_vecs(ik,3)*lz)*(0,1))*b_k**2
enddo
enddo


endif
enddo
enddo
enddo


if (ldiag) then	!on site energy; I ignore other on site elements

do i=0,3
Ham_k(i,i)=Ham_k(i,i)+Ham_input(0,0,0,i,i)-mu_k			!lambda and mu
enddo

do i=4,dim_hk-1
Ham_k(i,i)=Ham_k(i,i)+Ham_input(0,0,0,i,i)-mu_k-lambda_k
enddo

endif


	endif	!lreadham

!if (n_sites <5) then
do i=0, dim_hk-1
do j=i, dim_hk-1
#ifdef cmplx
if (abs(Ham_k(j,i)-conjg(Ham_k(i,j)))>1e-9) then
#else
if (abs(Ham_k(j,i)-(Ham_k(i,j)))>1e-9) then
#endif
print *, "error in Hamiltonian", ik, i, j, Ham_k(i,j), Ham_k(j,i)
stop
endif
enddo
enddo
!endif




end subroutine build_ham_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_param_kdotp
real(kind=idp) :: e_gamma, e_x, v_0, v_1, v_2,v_0_old, v_1_old,gl1d_t,gl1f_t,theta, theta_values(6),phi_values(6), ab_t, f1_t,v1_t,v2_t,E_t,eta3_f,B_0,B_1,B_t,v1_told,v2_told
 character(3):: label_theta(6)


if (phi_type=="smb" .or. phi_type=="sm7") then
f1=2*eta_v+2*etav2_x1+6*etav2_z2
f2=2*eta_vz+6*etav2_x1+2*etav2_z2
f1k2=-3*etav2_z2/2
h1k2=-3*etav2_z2/2
h1=-0.5*eta_v-1.5*eta_vz-3*etav2_x1+3*etav2_z2
h2=-1.5*eta_v-0.5*eta_vz+3*etav2_x1-3*etav2_z2
f7=2*eta_v7-4*sqrt(3.0d0)*etav2_z7-4*etav2_x7
h7=eta_v7+2*sqrt(3.0d0)*etav2_z7-6*etav2_x7
h72=sqrt(3.)*eta_v7-2*etav2_z7+2*sqrt(3.)*etav2_x7
hx=sqrt(3.0d0)/2*(-eta_v + eta_vz - 2*etav2_x1 + 2*etav2_z2)	!this must be split into hx12 and hx21 when including other hybridization terms
m78=(-4*eta_78 + 8*eta2_f781)*tf_*b_k**2
g1d=(eta1_d+3*eta2_d)
g2d=(eta+eta2_d)
g1f=(eta1_f+3*eta2_f82)
g2f=(eta_f+eta2_f82)
g7=(eta_7+4*eta2_f7+4*eta3_f7)
l7=(-eta_7+4*eta3_f7)
l1d=(-eta1_d-3*eta+6*eta2_d)/4
l2d=(-3*eta1_d-eta-6*eta2_d)/4
l1f=(-eta1_f-3*eta_f+6*eta2_f82)/4
l2f=(-3*eta1_f-eta_f-6*eta2_f82)/4
eckp1=ec_site_k(0)+tc_*(3*eta-eta1_d-6*eta2_d)
eckp2=ec_site_k(1)+tc_*(3*eta1_d-eta+2*eta2_d)
efkp1=ef_site_k(0)+tf_*b_k**2*(3*eta_f-eta1_f-6*eta2_f82)
efkp2=ef_site_k(1)+tf_*b_k**2*(3*eta1_f-eta_f+2*eta2_f82)
if (l_orbital_f==3) efkp7=ef_site_k(2)+tf_*b_k**2*(2*eta_7-4*eta2_f7-8*eta3_f7)
td12=sqrt(3.)/4*(eta1_d-eta-2*eta2_d)
tf12=sqrt(3.)/4*(eta1_f-eta_f-2*eta2_f82)

!f1=f1/2
!h1=h1/2
!f2=f2/2
!h2=h2/2
!hx=hx/2

if (m78>0) delta_c=m78*(f1**2-f7**2)/abs(f7*f1)
if (m78<0) delta_c=abs(m78)*(h1**2-h7**2)/abs(h7*h1)




elseif(phi_type=="sm2" .or. phi_type=="trv" .or. phi_type=="ons") then
eta1_f=eta_7
eta2_f=eta2_f7
eta3_f=eta3_f7

f1=-eta_v/2-8*eta_v2		!there is a /4 in the definition of the hopping parameter in the 1st nn
h1=eta_v/2
g1d=(eta1_d+4*eta2_d+4*eta3_d)	!kz^2
g1f=(eta1_f+4*eta2_f+4*eta3_f)	
l1d=(-eta1_d+4*eta3_d)		!kpar^2
l1f=(-eta1_f+4*eta3_f)
eckp1=ec_site_k(0)+tc_*(2*eta1_d-4*eta2_d-8*eta3_d)
efkp1=ef_site_k(0)+tf_*b_k**2*(2*eta1_f-4*eta2_f-8*eta3_f)
f1k2=2*eta_v2
h1k2=2*eta_v2
endif


!analytical results for gamma_8, and for 4band model
if (l_orbital_f<3) B_0=(eckp1-efkp1+lambda_k)/(tc_*g1d-tf_*b_k**2*g1f)
if (l_orbital_f<3) B_1=(eckp1-efkp1+lambda_k)/(tc_*l1d-tf_*b_k**2*l1f)
if (l_orbital_f<3) e_gamma=	((efkp1-mu_k)*tc_*g1d-(eckp1-mu_k)*tf_*b_k**2*g1f)/(tc_*g1d-tf_*b_k**2*g1f)
if (l_orbital_f<3) e_x=	((efkp1-mu_k)*tc_*l1d-(eckp1-mu_k)*tf_*b_k**2*l1f)/(tc_*l1d-tf_*b_k**2*l1f)
if (l_orbital_f<3) v_0_old=2*b_k*Vcf_*h1*sqrt(-tc_*tf_*b_k**2*g1d*g1f)/(tf_*b_k**2*g1f-tc_*g1d)
if (l_orbital_f<3) v_1_old=2*b_k*Vcf_*f1*sqrt(-tc_*tf_*b_k**2*l1d*l1f)/(tf_*b_k**2*l1f-tc_*l1d)
if (l_orbital_f<3) v_0=2*b_k*Vcf_*(h1-B_0*h1k2)*sqrt(-tc_*tf_*b_k**2*g1d*g1f)/(tf_*b_k**2*g1f-tc_*g1d)
if (l_orbital_f<3) v_1=2*b_k*Vcf_*(f1-B_1*f1k2)*sqrt(-tc_*tf_*b_k**2*l1d*l1f)/(tf_*b_k**2*l1f-tc_*l1d)
if (l_orbital_f<3) v_2=2*b_k*Vcf_*h1*sqrt(-tc_*tf_*b_k**2*l1d*l1f)/(tf_*b_k**2*l1f-tc_*l1d)
!analytical results for gamma_7
if (l_orbital_f==3) B_0=(eckp1-efkp7+lambda_k)/(tc_*g1d-tf_*b_k**2*g7)
if (l_orbital_f==3) B_1=(eckp1-efkp7+lambda_k)/(tc_*l1d-tf_*b_k**2*l7)
if (l_orbital_f==3) e_gamma=	((efkp7-mu_k)*tc_*g1d-(eckp1-mu_k)*tf_*b_k**2*g7)/(tc_*g1d-tf_*b_k**2*g7)
if (l_orbital_f==3) e_x=	((efkp7-mu_k)*tc_*l1d-(eckp1-mu_k)*tf_*b_k**2*l7)/(tc_*l1d-tf_*b_k**2*l7)
if (l_orbital_f==3) v_0=2*b_k*Vcf_*h7*sqrt(-tc_*tf_*b_k**2*g1d*g7)/(tf_*b_k**2*g7-tc_*g1d)
if (l_orbital_f==3) v_1=2*b_k*Vcf_*f7*sqrt(-tc_*tf_*b_k**2*l1d*l7)/(tf_*b_k**2*l7-tc_*l1d)
if (l_orbital_f==3) v_2=2*b_k*Vcf_*h7*sqrt(-tc_*tf_*b_k**2*l1d*l7)/(tf_*b_k**2*l7-tc_*l1d)




  OPEN(unit=725,file=trim(label)//'/par_kdotp',status='unknown')

write (725, "(a10,f15.7)"), "f1=", f1
write (725, "(a10,f15.7)"), "f2=", f2
write (725, "(a10,f15.7)"), "h1=", h1
write (725, "(a10,f15.7)"), "h2=", h2
write (725, "(a10,f15.7)"), "f7=", f7
write (725, "(a10,f15.7)"), "h7=", h7
write (725, "(a10,f15.7)"), "h72=", h72
write (725, "(a10,f15.7)"), "hx=", hx
write (725, "(a10,f15.7)"), "m78=", m78
write (725, "(a10,f15.7)"), "g1d=", g1d
write (725, "(a10,f15.7)"), "g2d=", g2d
write (725, "(a10,f15.7)"), "g1f=", g1f
write (725, "(a10,f15.7)"), "g2f=", g2f
write (725, "(a10,f15.7)"), "g7=", g7
write (725, "(a10,f15.7)"), "l7=", l7
write (725, "(a10,f15.7)"), "l1d=", l1d
write (725, "(a10,f15.7)"), "l2d=", l2d
write (725, "(a10,f15.7)"), "l1f=", l1f
write (725, "(a10,f15.7)"), "l2f=", l2f
write (725, "(a10,f15.7)"), "eckp1=", eckp1
write (725, "(a10,f15.7)"), "eckp2=", eckp2
write (725, "(a10,f15.7)"), "efkp1=", efkp1
write (725, "(a10,f15.7)"), "efkp2=", efkp2
write (725, "(a10,f15.7)"), "efkp7=", efkp7
write (725, "(a10,f15.7)"), "td12=", td12
write (725, "(a10,f15.7)"), "tf12=", tf12
write (725, "(a10,f15.7)"), "B_0=", B_0
write (725, "(a10,f15.7)"), "B_1=", B_1
write (725, "(a10,f15.7)"), "e7x-e1x=", efkp7-efkp1
write (725, "(a10,f15.7)"), "Delta_c=", delta_c
write (725, "(a10,f15.7)"), "e_gamma=", e_gamma
write (725, "(a10,f15.7)"), "e_x=", e_x
write (725, "(a10,f15.7)"), "v_0=", v_0
write (725, "(a10,f15.7)"), "v_1=", v_1
write (725, "(a10,f15.7)"), "v_2=", v_2
write (725, "(a10,f15.7)"), "v0old=", v_0_old
write (725, "(a10,f15.7)"), "v1old=", v_1_old


!other surfaces (just gamma8)
theta_values(1)=0	!001
theta_values(2)=pi/2	!100 110 210
theta_values(3)=pi/4	!101
theta_values(4)=atan(sqrt(2.))	!111
theta_values(5)=atan(2.)	!021
theta_values(6)=atan(1./2.)	!102

label_theta(1)="001"
label_theta(2)="100"
label_theta(3)="101"
label_theta(4)="111"
label_theta(5)="021"
label_theta(6)="102"

do i=1,6
theta=theta_values(i)

gl1d_t=g1d*cos(theta)**2+l1d*sin(theta)**2
gl1f_t=g1f*cos(theta)**2+l1f*sin(theta)**2
E_t= ((efkp1-mu_k)*tc_*gl1d_t-(eckp1-mu_k)*tf_*b_k**2*gl1f_t)/(tc_*gl1d_t-tf_*b_k**2*gl1f_t)
B_t=(eckp1-efkp1+lambda_k)/(tc_*gl1d_t-tf_*b_k**2*gl1f_t)


f1_t=sqrt(f1**2*cos(theta)**2+h1**2*sin(theta)**2)
!ab_t=sqrt(-tc_*tf_*b_k**2*gl1d_t*gl1f_t)/(tf_*b_k**2*gl1f_t-tc_*gl1d_t)*Vcf_*b_k*f1_t
v1_t=2*b_k*Vcf_* (h1*f1 -f1*h1k2*cos(theta)**2*(3*cos(theta)**2-2)*B_t-h1*f1k2*sin(theta)**2*(3*sin(theta)**2-2)*B_t)/f1_t  *sqrt(-tc_*tf_*b_k**2*gl1d_t*gl1f_t)/(tf_*b_k**2*gl1f_t-tc_*gl1d_t)
v2_t=2*b_k*Vcf_*(h1-B_t*cos(theta)**2*h1k2)	*sqrt(-tc_*tf_*b_k**2*gl1d_t*gl1f_t)/(tf_*b_k**2*gl1f_t-tc_*gl1d_t)
v1_told=2*b_k*Vcf_* (h1*f1)/f1_t  *sqrt(-tc_*tf_*b_k**2*gl1d_t*gl1f_t)/(tf_*b_k**2*gl1f_t-tc_*gl1d_t)
v2_told=2*b_k*Vcf_*(h1)	*sqrt(-tc_*tf_*b_k**2*gl1d_t*gl1f_t)/(tf_*b_k**2*gl1f_t-tc_*gl1d_t)


write (725, "(a3,a10,f15.7)"), ""
write (725, "(a3,a10,f15.7)"), label_theta(i),"theta=", theta_values(i)
write (725, "(a3,a10,f15.7)"), label_theta(i),"gl1d=", gl1d_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"gl1f=", gl1f_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"E_t=",E_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"B_t=",B_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"v1_t=",v1_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"v2_t=",v2_t
write (725, "(a3,a10,f15.7)"), label_theta(i),"v1t_old=",v1_told
write (725, "(a3,a10,f15.7)"), label_theta(i),"v2t_old=",v2_told
enddo


 close (unit=725)



end subroutine set_param_kdotp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kdotp(ik,Ham_k,b_k,lambda_k, mu_k,ldiag,nk,phi_sa, phi_as, ke, k_vecs)
complex(KIND=idpc),intent(out) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1)
complex(KIND=idpc),intent(in) 	:: phi_sa(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1),phi_as(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1),ke(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,0:1)
real(KIND=idp),intent(in) 	:: b_k, lambda_k, mu_k, k_vecs(0:nk-1,4)
integer, intent(in)		::ik,nk
logical,intent(in)		::ldiag
complex(KIND=idpc)		::kp,km
real(KIND=idp) 	:: D,Dxy,Dz,F1e,F1xy,F1z,F3,F3xy,F3z,F5,F5xy,F5z,c_1,c_2,c_3,c_4,c_5,d_1
complex(KIND=idpc),allocatable:: ham_k_135(:,:),rot_135_78(:,:)
!order of states: c up, c down, f +, f -

!eta_v=eta_x^v1
!eta_vz=eta_z^v1
!etav2_x1=eta_x^v2
!etav2_z2=eta_z^v2
!eta_v7=eta_7^v1
!etav2_z7=eta_7^v2
!etav2_x7=eta_x7^v2



kp=(k_vecs(ik,1)+(0,1)*k_vecs(ik,2))*Vcf_*b_k*2*pi	!kx+iky
km=(k_vecs(ik,1)-(0,1)*k_vecs(ik,2))*Vcf_*b_k*2*pi	!kx-iky
kz=(k_vecs(ik,3)-0.5)*Vcf_*b_k*2*pi			!kz^2
kpar2=(k_vecs(ik,1)**2+k_vecs(ik,2)**2)*4*pi**2		!kx^2+ky^2
kz2=(k_vecs(ik,3)-0.5)**2*4*pi**2			!(kz-pi)^2
k2xy=(-k_vecs(ik,1)**2+k_vecs(ik,2)**2)*4*pi**2		!ky^2-kx^2


if (ik==0) then
print *, "f1=",f1
print *, "f2=",f2
print *, "h1=",h1
print *, "h2=",h2
print *, "f7=",f7
print *, "h7=",h7
print *, "h72=",h72
print *, "hx=",hx
print *, "m78=",m78
endif


Ham_k=0

!eta=eta_z^d1
!eta1_d=eta_x^d1
!eta2_d=eta_z^d2
!eta_f=eta_z^f1
!eta1_f=eta_x^f1
!eta2_f82=eta_z^f2


if (.not. lreadham) then
!d1 d2 f1 f2 f7
!onsite energies
if (l_orbital==2) then
ham_k(0,0)=eckp1-mu_k+tc_*kz2*g1d+tc_*kpar2*l1d	!e_d1
ham_k(1,1)=ham_k(0,0)
ham_k(2,2)=eckp2-mu_k+tc_*kz2*g2d+tc_*kpar2*l2d	!e_d2
ham_k(3,3)=ham_k(2,2)
ham_k(4,4)=efkp1-mu_k-lambda_k+tf_*b_k**2*kz2*g1f+tf_*b_k**2*kpar2*l1f	!ef1
ham_k(5,5)=ham_k(4,4)
ham_k(6,6)=efkp2-mu_k-lambda_k+tf_*b_k**2*kz2*g2f+tf_*b_k**2*kpar2*l2f	!e_f2
ham_k(7,7)=ham_k(6,6)
if (l_orbital_f==3) ham_k(8,8)=efkp7-mu_k-lambda_k+tf_*b_k**2*kz2*g7+tf_*b_k**2*l7*kpar2	!e_f7
if (l_orbital_f==3) ham_k(9,9)=ham_k(8,8)
!non diag kin terms
ham_k(0,2)=tc_*td12*k2xy
ham_k(1,3)=ham_k(0,2)
ham_k(4,6)=tf_*b_k**2*tf12*k2xy
ham_k(5,7)=ham_k(4,6)

ham_k(0,4)=-(0,1)*f1*kz
ham_k(0,5)=-(0,1)*h1*km
ham_k(0,7)=-(0,1)*hx*kp
if (l_orbital_f==3) ham_k(0,8)=-(0,1)*f7*kz
if (l_orbital_f==3) ham_k(0,9)=-(0,1)*h7*km
ham_k(1,4)=-(0,1)*h1*kp
ham_k(1,5)=+(0,1)*f1*kz
ham_k(1,6)=-(0,1)*hx*km
if (l_orbital_f==3) ham_k(1,8)=-(0,1)*h7*kp
if (l_orbital_f==3) ham_k(1,9)=+(0,1)*f7*kz
ham_k(2,5)=-(0,1)*hx*kp
ham_k(2,6)=-(0,1)*f2*kz
ham_k(2,7)=-(0,1)*h2*km
if (l_orbital_f==3) ham_k(2,9)=-(0,1)*h72*kp
ham_k(3,4)=-(0,1)*hx*km
ham_k(3,6)=-(0,1)*h2*kp
ham_k(3,7)=+(0,1)*f2*kz
if (l_orbital_f==3) ham_k(3,8)=-(0,1)*h72*km
if (l_orbital_f==3) ham_k(4,8)=m78
if (l_orbital_f==3) ham_k(5,9)=m78

else	!4 band model
ham_k(0,0)=eckp1-mu_k+tc_*kz2*g1d+tc_*kpar2*l1d	!e_d1
ham_k(1,1)=ham_k(0,0)
ham_k(2,2)=efkp1-mu_k-lambda_k+tf_*b_k**2*kz2*g1f+tf_*b_k**2*kpar2*l1f	!ef1
ham_k(3,3)=ham_k(2,2)
ham_k(0,2)=-(0,1)*f1*kz
ham_k(1,3)=+(0,1)*f1*kz
ham_k(0,3)=-(0,1)*h1*km
ham_k(1,2)=-(0,1)*h1*kp
endif



do i=0,dim_hk-1
do j=i,dim_hk-1
ham_k(j,i)=conjg(ham_k(i,j))
enddo
enddo


endif

call write_Ham_k(ham_k)



!if (n_sites <5) then
do i=0, dim_hk-1
do j=i, dim_hk-1
#ifdef cmplx
if (abs(Ham_k(j,i)-conjg(Ham_k(i,j)))>1e-9) then
#else
if (abs(Ham_k(j,i)-(Ham_k(i,j)))>1e-9) then
#endif
print *, "error in Hamiltonian", ik, i, j, Ham_k(i,j), Ham_k(j,i)
stop
endif
enddo
enddo
!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!yu, xi dai
if (LREADHAM) then
kp=(k_vecs(ik,1)+(0,1)*k_vecs(ik,2))!*2*pi	!kx+iky
km=(k_vecs(ik,1)-(0,1)*k_vecs(ik,2))!*2*pi	!kx-iky
kz=(k_vecs(ik,3)-0.5)!*2*pi			!kz^2
kpar2=(k_vecs(ik,1)**2+k_vecs(ik,2)**2)!*4*pi**2		!kx^2+ky^2
kz2=(k_vecs(ik,3)-0.5)**2!*4*pi**2			!(kz-pi)^2


 c_1=-0.2787
 c_2=-0.318
 c_3=-0.502
 c_4=-0.5066
 c_5=0.0021
 d_1=-0.0132
! ed1=
! ef5=
! ef3=
! ef1=
 D=-1.5698
 Dxy=29.9233
 Dz=18.2502
 F1e=-0.0502
 F3=0.0193
 F5=-0.0134
 F1xy=-0.0020
 F3xy=-0.5776
 F5xy=-0.0532
 F1z=-0.0300
 F3z=-0.4476
 F5z=-0.2152
 
 
 allocate(ham_k_135(0:9,0:9))
 allocate(rot_135_78(0:9,0:9))
ham_k_135=0 
rot_135_78=0

rot_135_78(1,0)=1	!d_x^2-y2 down
rot_135_78(0,1)=1	!up
rot_135_78(3,8)=1	!d_z^2 down
rot_135_78(2,9)=1	!up
rot_135_78(5,3)=sqrt(5./6.)	!gamma_81 - 
rot_135_78(5,4)=sqrt(1./6.)
rot_135_78(4,2)=sqrt(5./6.)	!gamma_81+
rot_135_78(4,5)=sqrt(1./6.)
rot_135_78(7,7)=1	!gamma_82-
rot_135_78(6,6)=1	!gamma_82+
rot_135_78(9,3)=sqrt(1./6.)	!gamma_7 -
rot_135_78(9,4)=-sqrt(5./6.)	
rot_135_78(8,2)=sqrt(1./6.)	!gamma_7 +
rot_135_78(8,5)=-sqrt(5./6.)	

!onsite
ham_k_135(0,0)=D+Dxy*kpar2+Dz*kz2-2*delta_e7_read/6
ham_k_135(1,1)=ham_k_135(0,0)
ham_k_135(2,2)=F5+F5xy*kpar2+F5z*kz2-delta_e7_read/6
ham_k_135(3,3)=ham_k_135(2,2)
ham_k_135(4,4)=F3+F3xy*kpar2+F3z*kz2+3*delta_e7_read/6
ham_k_135(5,5)=ham_k_135(4,4)
ham_k_135(6,6)=F1e+F1xy*kpar2+F1z*kz2-2*delta_e7_read/6
ham_k_135(7,7)=ham_k_135(6,6)
ham_k_135(8,8)=50
ham_k_135(9,9)=ham_k_135(8,8)
!35
ham_k_135(3,4)=d_1-sqrt(5.)/6*+delta_e7_read
ham_k_135(4,3)=ham_k_135(3,4)
ham_k_135(5,2)=ham_k_135(3,4)
ham_k_135(2,5)=ham_k_135(3,4)
!hybr
ham_k_135(0,2)=c_1*kp
ham_k_135(1,3)=-c_1*km
ham_k_135(0,3)=c_2*kz
ham_k_135(1,2)=c_2*kz
ham_k_135(0,4)=c_3*kz
ham_k_135(1,5)=c_3*kz
ham_k_135(0,5)=c_4*kp
ham_k_135(1,4)=-c_4*km
ham_k_135(0,6)=c_5*km
ham_k_135(1,7)=-c_5*kp


ham_k_135(2,0)=conjg(ham_k_135(0,2))
ham_k_135(3,1)=conjg(ham_k_135(1,3))
ham_k_135(3,0)=conjg(ham_k_135(0,3))
ham_k_135(2,1)=conjg(ham_k_135(1,2))
ham_k_135(4,0)=conjg(ham_k_135(0,4))
ham_k_135(5,1)=conjg(ham_k_135(1,5))
ham_k_135(5,0)=conjg(ham_k_135(0,5))
ham_k_135(4,1)=conjg(ham_k_135(1,4))
ham_k_135(6,0)=conjg(ham_k_135(0,6))
ham_k_135(7,1)=conjg(ham_k_135(1,7))


print *, "ham Xi Dai:"
do i=0,9
write(*, "(20(f7.3))"), (ham_k_135(i,j),j=0,9)
enddo

!print *, "rotation:"
!do i=0,9
!write(*, "(20(f7.3))"), (rot_135_78(i,j),j=0,9)
!enddo


ham_k=0
ham_k=ham_k_135
!call BASIS_ROTATIONcmplx (transpose(rot_135_78), ham_k_135, transpose(rot_135_78), ham_k)

endif

end subroutine build_ham_kdotp




!!!!!!!!!!!!!!!!!!!!

subroutine build_ham_kdotp_red(Ham_k,Ham_k_red,rot_kdp)
complex(KIND=idpc),intent(in) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1),rot_kdp(0:dim_hk-1,0:3)
complex(KIND=idpc),intent(out) 	:: Ham_k_red(0:3,0:3)


!ham_k_red = MATMUL( rot_kdp, MATMUL(ham_k, transpose(rot_kdp)) )

call BASIS_ROTATIONcmplx_c (rot_kdp, ham_k, rot_kdp, ham_k_red,dim_hk,4)


end subroutine build_ham_kdotp_red

!!!!!!!!!!!!!!!!!!!!!
subroutine set_rot_kdotp(rot_kdp)
complex(KIND=idpc),intent(out) 	:: rot_kdp(0:dim_hk-1,0:3)
real(kind=idp)	::delta_17,tanxi

rot_kdp=0


rot_kdp(0,0)=1	!d1up
rot_kdp(1,1)=1	!d1down


if (dim_hk==4) then
rot_kdp(2,2)=1	!
rot_kdp(3,3)=1	!
elseif (dim_hk==6) then
rot_kdp(4,2)=1	!f7+
rot_kdp(5,3)=1	!f7-
elseif (dim_hk==8) then
rot_kdp(4,2)=1	!f1+
rot_kdp(5,3)=1	!f1-
elseif (dim_hk==10) then
delta_17=efkp7-efkp1
if (abs(m78)>1e-6) then
tanxi=delta_17/(2*m78)+m78/abs(m78)*sqrt((delta_17/(2*m78))**2+1)

!print *, efkp7, efkp1
!print *, "delta_17=", delta_17
!print *, "ep_m_e1=", ep_m_e1
print *, "m78", m78
print *, "tanxi", tanxi


!rot_kdp(4,2)=m78/sqrt(m78**2+ep_m_e1**2)	!f17+/f1+
!rot_kdp(8,2)=ep_m_e1/sqrt(ep_m_e1**2+m78**2)	!f17+/f7+
!rot_kdp(5,3)=m78/sqrt(m78**2+ep_m_e1**2)	!f17-/f1-
!rot_kdp(9,3)=ep_m_e1/sqrt(ep_m_e1**2+m78**2)	!f17-/f7-

rot_kdp(4,2)=-tanxi/sqrt(tanxi**2+1)	!f17+/f1+
rot_kdp(8,2)=1/sqrt(tanxi**2+1)	!f17+/f7+
rot_kdp(5,3)=-tanxi/sqrt(tanxi**2+1)	!f17-/f1-
rot_kdp(9,3)=1/sqrt(tanxi**2+1)	!f17-/f7-


else !(abs(m78)<1e-6) then
if (delta_17>0) then
rot_kdp(4,2)=1	!f17+/f1+
rot_kdp(5,3)=1	!f17-/f1-
else
rot_kdp(8,2)=1	!f17+/f7+
rot_kdp(9,3)=1	!f17-/f7-
endif
endif

endif

print *, "rot_kdp:"
do i=0,dim_hk-1
write(*, "(100(f10.5))"), (rot_kdp(i,j),j=0,3)
enddo





end subroutine set_rot_kdotp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_v_k(v_k_tot)
complex(kind=idp),intent(out):: v_k_tot(0:dim_hk-1,0:dim_hk-1, 0:nk3d-1)
integer	::ik

v_k_tot=0

do ik=0, nk3d-1
 call  build_ham_k(ik,V_k,1.d0,0.0d0, 0.0d0,.false.,nk3d,phi_k_sigma_alpha_tot,phi_k_alpha_sigma_tot, kin_energy_k_tot,k_vecs3d)
 v_k_tot(:,:,ik)=v_k(:,:)

!write (*, "(i5,3f10.5)") ik, k_vecs3d(ik,1),k_vecs3d(ik,2),k_vecs3d(ik,3)
!call write_ham_k(v_k)
enddo

end subroutine build_v_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sum_i(nn) exp(i*k*r_i)

subroutine factor_K(f_k_tot, n_kin, k_vecs,d23)
real (kind=idp),intent(out)	::f_k_tot(0:n_kin-1)
integer,intent(in)	::n_kin,d23
real (kind=idp),intent(in)	::k_vecs(0:n_kin-1,4)
integer			::i,ik,nn_max
complex(kind=idpc)	::f_k_c


if (d23==3) then
nn_max=n_nn
else if (d23==2) then
if (lkz) nn_max=n_nn-2
if (.not. lkz) nn_max=n_nn
else
print *, "error set_factor_k"
stop
endif



do ik=0, n_kin-1

f_k_c=0

do i=1, nn_max
!f_k_c=f_k_c+exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))
f_k_c=f_k_c+exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))


!print *, nn_ijk(i,1), k_vecs(ik,1), nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3),&
!(2*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))!,exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))

enddo


if (aimag(f_k_c)>1e-9) print *, "warning f_c_k", f_k_c

f_k_tot(ik)=real(f_k_c)

enddo

end subroutine factor_K


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_ham_k(Ham_k)
complex(KIND=idpc) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1)

do, i=0,dim_hk-1
    write(*,"(100(a1,f6.2,a1,f6.2,a1))") ( "(",real(Ham_k(i,j)),",",aimag(Ham_k(i,j)),")", j=0,dim_hk-1 )
enddo
   write(*, " ")


end subroutine write_ham_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_mat(Ham_k,dim_hk)
integer:: dim_hk
complex(KIND=idpc) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1)

do, i=0,dim_hk-1
    write(*,"(100(a1,f6.2,a1,f6.2,a1))") ( "(",real(Ham_k(i,j)),",",aimag(Ham_k(i,j)),")", j=0,dim_hk-1 )
enddo
   write(*, " ")


end subroutine write_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_mat2(Ham_k,dim_hk, dim_hk2)
integer:: dim_hk, dim_hk2
complex(KIND=idpc) 	:: Ham_k(0:dim_hk-1,0:dim_hk2-1)

do, i=0,dim_hk-1
    write(*,"(100(a1,f6.2,a1,f6.2,a1))") ( "(",real(Ham_k(i,j)),",",aimag(Ham_k(i,j)),")", j=0,dim_hk2-1 )
enddo
   write(*, " ")


end subroutine write_mat2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagonalize_ham_k(Ham_k, E_k, W_k)
complex(KIND=idpc) 	:: Ham_k(0:dim_hk-1,0:dim_hk-1),W_k(0:dim_hk-1,0:dim_hk-1)
real(KIND=idp)		:: E_k(0:dim_hk-1)

CALL DIAGONALIZE( 'V', dim_hk, dim_hk, Ham_k, E_k, W_k )


end subroutine diagonalize_ham_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_energies_k(e_k,w_k)
real(KIND=idp) 	:: E_k(0:dim_hk-1)
complex(KIND=idpc) 	:: w_k(0:dim_hk-1,0:dim_hk-1)

do, i=0,dim_hk-1
    write(*,"(100f15.10)") ( e_k(i) )
enddo

!do, i=0,3
!   write(*,"(100f10.5)") ( real(w_k(i,j)),aimag(w_k(i,j)), j=0,3 )
!enddo



END SUBROUTINE write_energies_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine compute_mu_k(e_k_tot)
real(KIND=idp), intent(in)	:: E_k_tot(0:dim_hk-1, 0:nk3d-1)
integer			        :: i,j,ik
real (kind=idp)			::thresh, alfa_mu, dndmuc


!mu_old_k=mu_k

!i_iter_mu=0
!diff_mu=100
 
 
!write(*,"(a40,100f15.10)") "******* Mu loop, starting mu=", mu_k
!write(*,"(a5,5(a20),a20,a15)") "i_mu","mu",  "Nfctot",  "N_el", "diff_mu", "mu_new", "dndmu"

!!!!!!!mu loop


!if (i_iter<5) then
!thresh=0.1
!alfa_mu=0.005
!else
!thresh=thr_mu
!alfa_mu=alpha_mu_k
!endif

!alfa_mu=0.01*(1-1/(i_iter+0.1))
!thresh=thr_mu+0.01/i_iter

!alfa_mu=alpha_mu_k
!thresh=thr_mu



!do while (diff_mu>thresh)
!i_iter_mu=i_iter_mu+1

   
   
!   print *, e_k_tot

dndmu=0
nfctot=0
do i=0, dim_hk-1	!eigenvalues
 do ik=0,nk3d-1
!  nfctot=nfctot+1/(exp((e_k_tot(i,ik)-mu_k)/T)+1)*k_vecs3d(ik,4)
!  dndmu=dndmu+1/T/(exp((e_k_tot(i,ik)-mu_k)/T)+1)/(exp((-e_k_tot(i,ik)+mu_k)/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
  nfctot=nfctot+1/(exp((e_k_tot(i,ik))/T)+1)*k_vecs3d(ik,4)
  dndmu=dndmu+1/T/(exp((e_k_tot(i,ik))/T)+1)/(exp((-e_k_tot(i,ik))/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
 enddo
enddo

!dndmuc=0
!nctot=0

!      do ik=0,nk3d-1	!k -points
!       do i=0,4*l_orbital-1		!eigenvalues
!     nctot= nctot + ((abs(w_k_tot(0,i,ik)))**2+abs(w_k_tot(1,i,ik))**2) / (exp((e_k_tot(i,ik)-mu_k)/T)+1) * k_vecs3d(ik,4)
!     dndmuc=dndmuc+((abs(w_k_tot(0,i,ik)))**2+abs(w_k_tot(1,i,ik))**2)/T/(exp((e_k_tot(i,ik)-mu_k)/T)+1)/(exp((-e_k_tot(i,ik)+mu_k)/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
!     enddo
!    enddo



!print *, "dndmu",i_iter,dndmu,nfctot,n_el_k
dndmu=max(dndmu,0.001)		!to avoid dndmu --> 0 when T-->0
!dndmuc=max(dndmuc,0.001)		!to avoid dndmu --> 0 when T-->0

!!!!!!Difference between wanted and obtained N_el

!diff_mu=abs(real(N_el_k)-Nfctot)

!!!!!!!!!!!!! Change mu



!print *, dndmu, nfctot, n_el_k,(N_el_k-Nfctot)/dndmu


!if (abs(N_el_k-Nfctot)<N_k*0.1)then
!else
!dndmu=-N_k/(4*tc_)
!endif

if (lmucf) mu_new_k=mu_k+(N_el_k-Nfctot)/dndmu
!print *, dndmu, mu_new_k,lmucf
!print *, mu_k,mu_new_k
!if (.not. lmucf) mu_new_k=mu_k+(N_el_k/2-Nctot)/dndmuc

!print *, Nfctot,mu_k,mu_new_k

!print *, dndmu, mu_new_k

!if (mod(i_iter_mu,1001)==0) write(*,"(i6,5f20.10,f20.10, f15.10)") i_iter_mu, mu_k,  Nfctot,  real(N_el_k), diff_mu, mu_new_k, dndmu

!mu_new_k=alfa_mu*mu_new_k+(1-alfa_mu)*mu_k
!mu_k=mu_new_k

!enddo


!!!!!!!!!!Finish loop

!mu_k=mu_old_k


!write(*,"(a60,f15.10, i10, f15.5)") "********Mu loop finished, final mu, n_iter, n_fctot=", mu_new_k, i_iter_mu, Nfctot

!pause
end subroutine compute_mu_k



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_mu_k_tot(e_k_tot,mu_old_k)
real(KIND=idp), intent(in)	:: E_k_tot(0:dim_hk-1, 0:nk3d-1),mu_old_k
integer			        :: i,j,ik
real (kind=idp)			::thresh, alfa_mu, dndmuc,weight,tot_weight,n_tot

mu_k=mu_old_k
!mu_old_k=mu_k

!i_iter_mu=0
!diff_mu=100
 
 
write(*,"(a40,100f15.10)") "******* Mu loop, starting mu=", mu_k
write(*,"(a5,5(a20),a20,a15)") "i_mu","mu",  "Nfctot",  "N_el", "diff_mu", "mu_new", "dndmu"

!!!!!!!mu loop



!alfa_mu=0.01*(1-1/(i_iter+0.1))
!thresh=thr_mu+0.01/i_iter

!alfa_mu=alpha_mu_k
!thresh=thr_mu


diff_mu=10
i_iter_mu=0
thresh=thr_mu

do while (diff_mu>thresh)
i_iter_mu=i_iter_mu+1


if (i_iter<100) then
alfa_mu=alpha_mu_k*0.001
else
alfa_mu=alpha_mu_k
endif

   
   
!   print *, e_k_tot

dndmu=0
nfctot=0
do i=0, dim_hk-1	!eigenvalues
 do ik=0,nk3d-1
!  nfctot=nfctot+1/(exp((e_k_tot(i,ik)-mu_k)/T)+1)*k_vecs3d(ik,4)
!  dndmu=dndmu+1/T/(exp((e_k_tot(i,ik)-mu_k)/T)+1)/(exp((-e_k_tot(i,ik)+mu_k)/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
  nfctot=nfctot+1/(exp((e_k_tot(i,ik)-mu_k)/T)+1)*k_vecs3d(ik,4)
  dndmu=dndmu+1/T/(exp((e_k_tot(i,ik)-mu_k)/T)+1)/(exp((-e_k_tot(i,ik)+mu_k)/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
 enddo
enddo
!print *, dndmu

!dndmuc=0
!nctot=0

!      do ik=0,nk3d-1	!k -points
!       do i=0,4*l_orbital-1		!eigenvalues
!     nctot= nctot + ((abs(w_k_tot(0,i,ik)))**2+abs(w_k_tot(1,i,ik))**2) / (exp((e_k_tot(i,ik)-mu_k)/T)+1) * k_vecs3d(ik,4)
!     dndmuc=dndmuc+((abs(w_k_tot(0,i,ik)))**2+abs(w_k_tot(1,i,ik))**2)/T/(exp((e_k_tot(i,ik)-mu_k)/T)+1)/(exp((-e_k_tot(i,ik)+mu_k)/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
!     enddo
!    enddo



!print *, "dndmu",i_iter,dndmu,nfctot,n_el_k
dndmu=max(dndmu,0.001)		!to avoid dndmu --> 0 when T-->0
!dndmuc=max(dndmuc,0.001)		!to avoid dndmu --> 0 when T-->0

!!!!!!Difference between wanted and obtained N_el

diff_mu=abs(real(N_el_k)-Nfctot)

!!!!!!!!!!!!! Change mu



!print *, dndmu, nfctot, n_el_k,(N_el_k-Nfctot)/dndmu


!if (abs(N_el_k-Nfctot)<N_k*0.1)then
!else
!dndmu=-N_k/(4*tc_)
!endif

mu_new_k=mu_k+(N_el_k-Nfctot)/dndmu

!print *, mu_k,mu_new_k
!if (.not. lmucf) mu_new_k=mu_k+(N_el_k/2-Nctot)/dndmuc

!print *, Nfctot,mu_k,mu_new_k

!print *, dndmu, mu_new_k

!write(*,"(i6,5f20.10,f20.10, f15.10)") i_iter_mu, mu_k,  Nfctot,  real(N_el_k), diff_mu, mu_new_k, dndmu
if (mod(i_iter_mu,1001)==0 .or. i_iter_mu<10) write(*,"(i7,5f20.10,f20.10, f15.10)") i_iter_mu, mu_k,  Nfctot,  real(N_el_k), diff_mu, mu_new_k, dndmu

mu_new_k=alfa_mu*mu_new_k+(1-alfa_mu)*mu_k
mu_k=mu_new_k

enddo


!!!!!!!!!!Finish loop

!mu_k=mu_old_k
tot_weight=0
nc_k=0
nf_k=0
n_tot=0

      do ik=0,nk3d-1	!k -points
    tot_weight=tot_weight+ k_vecs3d(ik,4)
           do i=0,dim_hk-1		!eigenvalues
       
    weight=1/(exp((e_k_tot(i,ik)-mu_k)/T)+1) * k_vecs3d(ik,4)			
			
  
 !n_c
       do j=0,2*l_orbital-1				!0,1,2,3
        nc_k= nc_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo
 !n_f      
       do j=2*l_orbital,dim_hk-1			!4,5,6,7(,8,9)
        nf_k= nf_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo

       do j=2*l_orbital,dim_hk-1			!4,5,6,7(,8,9)
        n_tot= n_tot +  weight
       enddo


enddo
enddo

write(*,"(a60,f15.10, i10, f15.5)") "********Mu loop finished, final mu, n_iter, n_fctot=", mu_new_k, i_iter_mu, Nfctot
  OPEN(unit=131,file=trim(label)//'/nf',status='unknown')
write(131,"(a10,2f20.10)"), "tot_w=",tot_weight
write(131,"(a10,2f20.10)"), "n_c=",nc_k,nc_k/tot_weight
write(131,"(a10,2f20.10)"),"n_f=",nf_k,nf_k/tot_weight
write(131,"(a10,2f20.10)"),"n_c+n_f=",nc_k+nf_k,(nc_k+nf_k)/tot_weight
write(131,"(a10,2f20.10)"),"n_tot=",n_tot,(n_tot)/tot_weight
write(131,"(a10,2f20.10)"),"b=",sqrt(1-nf_k/tot_weight)
 close(unit=131)

!pause
end subroutine compute_mu_k_tot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_nf_k (nf_k,nc_k,b_new_k,e_k_tot,w_k_tot,ef_k, lambda_new_k, nfc_k, nff_k,diff, mu)
complex(KIND=idpc),intent(in) 	:: w_k_tot(0:dim_hk-1,0:dim_hk-1, 0:nk3d-1)
real(KIND=idp),intent(out) 		:: b_new_k,lambda_new_k, diff
real(KIND=idp),intent(in) 		:: E_k_tot(0:dim_hk-1, 0:nk3d-1), mu
real(kind=idp)		:: nf_k,nc_k,nfc_k,nff_k
real(kind=idp)		:: ec_k,ef_k, weight
integer		::i,j,ik,is,is2

nf_k=0
nf7_k=0
nc_k=0
nfc_k=0	!V included
nff_k=0	!tf included

  
      do ik=0,nk3d-1	!k -points
       do i=0,dim_hk-1		!eigenvalues
       
    weight=1/(exp((e_k_tot(i,ik)-mu)/T)+1) * k_vecs3d(ik,4)			!mu=0 (already in Hamiltonian)
  
 !n_c
       do j=0,2*l_orbital-1				!0,1,2,3
        nc_k= nc_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo
 !n_f      
       do j=2*l_orbital,dim_hk-1			!4,5,6,7(,8,9)
        nf_k= nf_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo
if (l_orbital_f==3) then
       do j=8,9			
        nf7_k= nf7_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo
endif
 !<f+f> including tf_
       do j=2*l_orbital,dim_hk-1
       do l=2*l_orbital,dim_hk-1
        nff_k= nff_k + conjg(w_k_tot(j,i,ik))*V_k_tot(j,l,ik)*w_k_tot(l,i,ik)* weight
       enddo
       enddo


!     nff_k=nff_k+ factor_k_tot(ik)*tf_*((abs(w_k_tot(2,i,ik)))**2+abs(w_k_tot(3,i,ik))**2) * weight


 !<f+c> including V

       do j=0,2*l_orbital-1			!d
       do l=2*l_orbital,dim_hk-1		!f
        nfc_k= nfc_k + conjg(w_k_tot(j,i,ik))*V_k_tot(j,l,ik)*w_k_tot(l,i,ik)* weight
        nfc_k= nfc_k + conjg(w_k_tot(l,i,ik))*V_k_tot(l,j,ik)*w_k_tot(j,i,ik)* weight
       enddo
       enddo


!        do is=0,1		!sigma
!         do is2=0,1	!alpha


!       nfc_k =nfc_k+real ( conjg(w_k_tot(is,i,ik))*V_k_tot(is,2+is2,ik)*w_k_tot(2+is2,i,ik)&
!                           +conjg(w_k_tot(2+is2,i,ik))*(V_k_tot(2+is2,is,ik))*w_k_tot(is,i,ik) ) * weight
!                           +conjg(w_k_tot(2+is2,i,ik))*conjg(V_k_tot(is,2+is2,ik))*w_k_tot(is,i,ik) ) * weight
	
!        enddo
!       enddo

 
     enddo
    enddo


nf_k=nf_k/n_k		!divide by ALL k points!
nf7_k=nf7_k/n_k		!divide by ALL k points!
nc_k=nc_k/n_k
nff_k=nff_k/n_k
nfc_k=nfc_k/n_k

!!!!!!!!!total free energy
E_tot_k=0

!from eigenvalues
do i=0,dim_hk-1
 do j=0,nk3d-1
  E_tot_k=E_tot_k+e_k_tot(i,j)/(exp((e_k_tot(i,j)-mu)/T)+1) * k_vecs3d(j,4)
 enddo
enddo

!from b
E_tot_k=E_tot_k-n_k*lambda_new_k*(b_new_k**2-1)

!from mu
E_tot_k=E_tot_k+mu*N_el_k

!from lambda**2
E_tot_k_l=E_tot_k-n_k*lambda_new_k**2/(2*lambda_par_k)

E_tot_k=E_tot_k/n_k
E_tot_k_l=E_tot_k_l/n_k

!!!!!!


!!!!!!set new variables
!SB
if(scf_type=='sb') then

b_new_k=(nfc_k/(lambda_k-nff_k)/2)
lambda_new_k=lambda_k+lambda_par_k*(1-b_new_k**2-nf_k)

 ef_k=ef_-lambda_new_k!-mu

!print *, nctot, N_el_k/2

diff=0
if ( .not. lfixmu_k .and. lmucf)  	diff=(1-b_new_k**2-nf_k)**2+((Nfctot-N_el_k))**2+ (b_new_k-nfc_k/(lambda_new_k-nff_k)/2)**2
!if ( .not. lfixmu_k .and. .not. lmucf) 	diff=(1-b_new_k**2-nf_k)**2+((Nctot-N_el_k/2))**2+(b_new_k-nfc_k/(lambda_new_k-nff_k)/2)**2
if (lfixmu_k)  		diff=(1-b_new_k**2-nf_k)**2+(b_new_k-(nfc_k/(lambda_new_k-nff_k)/2))**2
diff=sqrt(diff)

!GA
!in this case b=R=sqrt((1-P)/(1-P/2))
else if(scf_type=='ga') then

b_new_k=nfc_k/2/(-nff_k+2*lambda_k/(2-b_k**2)**2)
lambda_new_k=lambda_k+lambda_par_k*((1-b_new_k**2)/(1-b_new_k**2/2)-nf_k)

!if (nf_k<1) b_new_k=sqrt((1-nf_k)/(1-nf_k/2))
!if (nf_k>=1) b_new_k=1e-3
!lambda_new_k=(2*b_new_k*nff_k+nfc_k)/dpdr(b_new_k)

 ef_k=ef_-lambda_new_k!-mu

!print *, nctot, N_el_k/2

diff=0
if ( .not. lfixmu_k .and. lmucf)  	diff=(p_r(b_new_k)-nf_k)**2+((Nfctot-N_el_k))**2+ (2*nff_k*b_new_k+lambda_new_k*dpdr(b_k)+nfc_k)**2
print *, (p_r(b_new_k)-nf_k)**2, ((Nfctot-N_el_k))**2, (2*nff_k*b_new_k+lambda_new_k*dpdr(b_k)+nfc_k)**2
!if ( .not. lfixmu_k .and. .not. lmucf) 	diff=(1-b_new_k**2-nf_k)**2+((Nctot-N_el_k/2))**2+(b_new_k-nfc_k/(lambda_new_k-nff_k)/2)**2
!if (lfixmu_k)  		diff=(1-b_new_k**2-nf_k)**2+(b_new_k-(nfc_k/(lambda_new_k-nff_k)/2))**2
diff=sqrt(diff)

endif


!!!!!!!!!!Writes new parameters
if (i_iter==1 .or. (mod(i_iter,30)==0)) &
write (*, "(a5,a16,21(a12))"), "i","diff", "b_old_k","nf_k","b_new_k","bnfc","ef_k","bk^2+nfk","nfc_k","lambda","l_new", "mu_newk", "nc_k", "nf_k+nc_k", "nff_k", "nf7_k"
write(*,"(i6,f17.10,100f12.5)")  i_iter,diff, b_k, nf_k, b_new_k,(nfc_k/(lambda_k-nff_k)/2), ef_k, b_new_k**2+nf_k, nfc_k, lambda_k, lambda_new_k, mu_new_k, nc_k, nf_k+nc_k, nff_k, nf7_k


  OPEN(unit=631,file=trim(label)//'/b_k3',status='unknown')
  write(631,"(3f20.10)")  b_k, lambda_k,  mu_k	!only if I did a scf run
  close(unit=631)

end subroutine compute_nf_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_nfnc_k (nf_k,nc_k,b_new_k,e_k_tot,w_k_tot,ef_k, lambda_new_k, nfc_k, nff_k,diff, mu)
complex(KIND=idpc),intent(in) 	:: w_k_tot(0:dim_hk-1,0:dim_hk-1, 0:nk3d-1)
real(KIND=idp),intent(out) 		:: b_new_k,lambda_new_k, diff
real(KIND=idp),intent(in) 		:: E_k_tot(0:dim_hk-1, 0:nk3d-1), mu
real(kind=idp)		:: nf_k,nc_k,nfc_k,nff_k
real(kind=idp)		:: ec_k,ef_k, weight,dndmuf
integer		::i,j,ik,is,is2

nf_k=0
nf7_k=0
nc_k=0

dndmuf=0  
dndmu=0
nfctot=0
      do ik=0,nk3d-1	!k -points
       do i=0,dim_hk-1		!eigenvalues
       
    weight=1/(exp((e_k_tot(i,ik)-mu)/T)+1) * k_vecs3d(ik,4)			!mu=0 (already in Hamiltonian)

!  nfctot=nfctot+1/(exp((e_k_tot(i,ik))/T)+1)*k_vecs3d(ik,4)
!  dndmu=dndmu+1/T/(exp((e_k_tot(i,ik))/T)+1)/(exp((-e_k_tot(i,ik))/T)+1)*k_vecs3d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
  
 !n_c
       do j=0,2*l_orbital-1				!0,1,2,3
        nc_k= nc_k + abs(w_k_tot(j,i,ik))**2 * weight
       enddo
 !n_f      
       do j=2*l_orbital,dim_hk-1			!4,5,6,7(,8,9)
        nf_k= nf_k + abs(w_k_tot(j,i,ik))**2 * weight
        dndmuf=dndmuf+1/T/(exp((e_k_tot(i,ik))/T)+1)/(exp((-e_k_tot(i,ik))/T)+1)*k_vecs3d(ik,4)*abs(w_k_tot(j,i,ik))**2		!d(n)/d(mu)=sum_i f'(e_i-mu)
       enddo
!if (l_orbital_f==3) then
!       do j=8,9			
!        nf7_k= nf7_k + abs(w_k_tot(j,i,ik))**2 * weight
!       enddo
!endif
 

     enddo
    enddo


dndmuf=max(dndmuf,0.001)		!to avoid dndmu --> 0 when T-->0
dndmu=max(dndmu,0.001)		!to avoid dndmu --> 0 when T-->0
nfctot=nfctot/n_k		!divide by ALL k points!
nf_k=nf_k/n_k		!divide by ALL k points!
nc_k=nc_k/n_k
!nf7_k=nf7_k/n_k		!divide by ALL k points!


!!!!!!set new variables
b_new_k=b_k
lambda_new_k=lambda_k+(nf_0-nf_k)/dndmuf
!mu_new_k=mu_k+(nf_0+nc_0-nf_k-nc_k)/dndmu

diff=abs(nf_0-nf_k)+abs(nc_0-nc_k)

!!!!!!!!!!Writes new parameters
if (i_iter==1 .or. (mod(i_iter,10000)==0)) write (*, "(a5,a16,21(a10))"), "i","diff", "nf_k","nf_0","nc_k","nc_0","nfctot", "nf_0+nc_0", "lambda","mu"
if (i_iter<20 .or.  (mod(i_iter,1000)==0))  write(*,"(i6,f17.10,100f10.5)")  i_iter,diff, nf_k,nf_0,nc_k,nc_0, nfctot, nf_0+nc_0,lambda_new_k,mu_k


  OPEN(unit=631,file=trim(label)//'/b_k3',status='unknown')
  write(631,"(3f20.10)")  b_k, lambda_k,  mu_k	!only if I did a scf run
  close(unit=631)

end subroutine compute_nfnc_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_bands_k(e_k_tot)
real (kind=idp)		:: e_k_tot(0:dim_hk-1,0:nk3d-1)
!real(kind=idp)intent(in)		::k_vecs(0:n_k-1,1:3)
!real (kind=idp), allocatable		:: dos(:,:)
complex(KIND=idpc) 	:: w_k_tot(0:dim_hk-1,0:dim_hk-1,0:nk3d-1)
integer		::i,j,ik



  write(14,"(3a10, 5a15)") "#k_x","k_y","k_z", "e1", "e2", "e3", "e4"!, "Delta_k^2"

do ik= 0, nk3d-1

 write(14,"(3f10.5, 5f15.10)") k_vecs3d(ik,1),k_vecs3d(ik,2),k_vecs3d(ik,3), e_k_tot(0,ik), e_k_tot(1,ik), e_k_tot(2,ik), e_k_tot(3,ik)!, delta_k

enddo    



end subroutine write_bands_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine write_k_points(k_vecs3d)
real(kind=idp), intent(in)		::k_vecs3d(0:nk3d-1,4)
integer		::i,ik

!do i=6,24,18

 write(*,"(a15)") ""
 write(*,"(a20)")  "List of K-points:"
 write(*,"(a6,6a10)")  "ind", "kx", "ky", "kz", "weight", "F_k", "Delta_k^2"
 write(*,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
 write(*,"(a15)") ""


 do ik=0, nk3d-1 

  write(*,"(i8,6f10.5)")  ik, 	 k_vecs3d(ik,1), k_vecs3d(ik,2), k_vecs3d(ik,3),k_vecs3d(ik,4), delta_k_tot(ik)


if (l_orbital==2 .and. l_orbital_f==2) then

   write(*,"(a35)") "Kin_energy_d"

    do o1=0,1
    do is=0,1
 write(*,"(6 ('(',f6.2,',',f6.2,')','  ') )")  	real(kin_energy_k_tot(o1,0,is,0,ik,0)),aimag(kin_energy_k_tot(o1,0,is,0,ik,0)), &
 						real(kin_energy_k_tot(o1,0,is,1,ik,0)),aimag(kin_energy_k_tot(o1,0,is,1,ik,0)),&
   						real(kin_energy_k_tot(o1,1,is,0,ik,0)),aimag(kin_energy_k_tot(o1,1,is,0,ik,0)), &
   						real(kin_energy_k_tot(o1,1,is,1,ik,0)),aimag(kin_energy_k_tot(o1,1,is,1,ik,0))
    enddo 
    enddo


   write(*,"(a35)") "Kin_energy_f"

    do o1=0,1
    do is=0,1
 write(*,"(6 ('(',f6.2,',',f6.2,')','  ') )")  	real(kin_energy_k_tot(o1,0,is,0,ik,1)),aimag(kin_energy_k_tot(o1,0,is,0,ik,1)), &
 						real(kin_energy_k_tot(o1,0,is,1,ik,1)),aimag(kin_energy_k_tot(o1,0,is,1,ik,1)),&
   						real(kin_energy_k_tot(o1,1,is,0,ik,1)),aimag(kin_energy_k_tot(o1,1,is,0,ik,1)), &
   						real(kin_energy_k_tot(o1,1,is,1,ik,1)),aimag(kin_energy_k_tot(o1,1,is,1,ik,1))
    enddo 
    enddo




   write(*,"(a35)") "Phi_k_sigma_alpha"

    do o1=0,1
    do is=0,1
 write(*,"(6 ('(',f6.2,',',f6.2,')','  ') )")  	real(phi_k_sigma_alpha_tot(o1,0,is,0,ik)),aimag(phi_k_sigma_alpha_tot(o1,0,is,0,ik)), &
 						real(phi_k_sigma_alpha_tot(o1,0,is,1,ik)),aimag(phi_k_sigma_alpha_tot(o1,0,is,1,ik)),&
   						real(phi_k_sigma_alpha_tot(o1,1,is,0,ik)),aimag(phi_k_sigma_alpha_tot(o1,1,is,0,ik)), &
   						real(phi_k_sigma_alpha_tot(o1,1,is,1,ik)),aimag(phi_k_sigma_alpha_tot(o1,1,is,1,ik))
    enddo 
    enddo


   write(*,"(a35)") "Phi_k_alpha_sigma"

    do o1=0,1
    do is=0,1
 write(*,"(6 ('(',f6.2,',',f6.2,')','  ') )")  	real(phi_k_alpha_sigma_tot(o1,0,is,0,ik)),aimag(phi_k_alpha_sigma_tot(o1,0,is,0,ik)), &
 						real(phi_k_alpha_sigma_tot(o1,0,is,1,ik)),aimag(phi_k_alpha_sigma_tot(o1,0,is,1,ik)),&
   						real(phi_k_alpha_sigma_tot(o1,1,is,0,ik)),aimag(phi_k_alpha_sigma_tot(o1,1,is,0,ik)), &
   						real(phi_k_alpha_sigma_tot(o1,1,is,1,ik)),aimag(phi_k_alpha_sigma_tot(o1,1,is,1,ik))
    enddo 
    enddo




   write(*,"(a35)") "All"

    do o1=0,1
    do is=0,1
 write(*,"(8 ('(',f6.2,',',f6.2,')','  ') )")  	real(kin_energy_k_tot(o1,0,is,0,ik,0)),aimag(kin_energy_k_tot(o1,0,is,0,ik,0)), &
 						real(kin_energy_k_tot(o1,0,is,1,ik,0)),aimag(kin_energy_k_tot(o1,0,is,1,ik,0)),& 	
 						real(kin_energy_k_tot(o1,1,is,0,ik,0)),aimag(kin_energy_k_tot(o1,1,is,0,ik,0)), &
 						real(kin_energy_k_tot(o1,1,is,1,ik,0)),aimag(kin_energy_k_tot(o1,1,is,1,ik,0)),&
 						real(phi_k_sigma_alpha_tot(o1,0,is,0,ik)),aimag(phi_k_sigma_alpha_tot(o1,0,is,0,ik)), &
 						real(phi_k_sigma_alpha_tot(o1,0,is,1,ik)),aimag(phi_k_sigma_alpha_tot(o1,0,is,1,ik)),&
 						real(phi_k_sigma_alpha_tot(o1,1,is,0,ik)),aimag(phi_k_sigma_alpha_tot(o1,1,is,0,ik)), &
 						real(phi_k_sigma_alpha_tot(o1,1,is,1,ik)),aimag(phi_k_sigma_alpha_tot(o1,1,is,1,ik))
    enddo 
    enddo
    do o1=0,1
    do is=0,1
 write(*,"(8 ('(',f6.2,',',f6.2,')','  ') )")  	real(phi_k_alpha_sigma_tot(o1,0,is,0,ik)),aimag(phi_k_alpha_sigma_tot(o1,0,is,0,ik)), &
 						real(phi_k_alpha_sigma_tot(o1,0,is,1,ik)),aimag(phi_k_alpha_sigma_tot(o1,0,is,1,ik)),&
 						real(phi_k_alpha_sigma_tot(o1,1,is,0,ik)),aimag(phi_k_alpha_sigma_tot(o1,1,is,0,ik)), &
 						real(phi_k_alpha_sigma_tot(o1,1,is,1,ik)),aimag(phi_k_alpha_sigma_tot(o1,1,is,1,ik)),&
						real(kin_energy_k_tot(o1,0,is,0,ik,1)),aimag(kin_energy_k_tot(o1,0,is,0,ik,1)), &
 						real(kin_energy_k_tot(o1,0,is,1,ik,1)),aimag(kin_energy_k_tot(o1,0,is,1,ik,1)),& 	
 						real(kin_energy_k_tot(o1,1,is,0,ik,1)),aimag(kin_energy_k_tot(o1,1,is,0,ik,1)), &
 						real(kin_energy_k_tot(o1,1,is,1,ik,1)),aimag(kin_energy_k_tot(o1,1,is,1,ik,1))
    enddo 
    enddo


    						 
    write(*, " ")




endif



  write(*, " ")
 enddo

 write(*,"(a10)") "**********"

!enddo

end subroutine write_k_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_nn_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)
!real(kind=idp)		:: nn_theta_phi(N_nn,2),nn_xyz(N_nn,3)
!integer, intent(in)	:: nn_ijk(N_nn,3)
!complex(kind=idpc)	::phi_sigma_alpha(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn), phi_alpha_sigma(0:l_orbital-1,0:l_orbital-1,0:1,0:1, n_nn)
!complex(kind=idpc)	::kin_energy(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn,0:1)
integer	:: i,j


if (l_orbital==1 .or. l_orbital==0) then
do j=6,24,18

    write(j,"(a30,i5)") "Number of n.n.:", n_nn
    write(j,"(a60)") "Nearest neighbours (c wrt to f, which is at the origin):"
    write(j,"(3a6, 5a10)") "i", "j", "k", "theta", "phi", "x", "y", "z"
    write(j,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
    write(j,"(a15)") ""

do, i=1,n_nn
    write(j,"(3i6, 5f10.5)")  nn_ijk(i,1),nn_ijk(i,2),nn_ijk(i,3),nn_theta_phi(i,1),nn_theta_phi(i,2),nn_xyz(i,1),nn_xyz(i,2),nn_xyz(i,3)

if (l_orbital==1) then
print *, "phi_sa"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_sigma_alpha(0,0,is,0,i)),aimag(phi_sigma_alpha(0,0,is,0,i)), real(phi_sigma_alpha(0,0,is,1,i)),aimag(phi_sigma_alpha(0,0,is,1,i))
    enddo

print *, "phi_as"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_alpha_sigma(0,0,is,0,i)),aimag(phi_alpha_sigma(0,0,is,0,i)), real(phi_alpha_sigma(0,0,is,1,i)),aimag(phi_alpha_sigma(0,0,is,1,i))
    enddo

print *, "kin_s"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy(0,0,is,0,i,0)),aimag(kin_energy(0,0,is,0,i,0)), real(kin_energy(0,0,is,1,i,0)),aimag(kin_energy(0,0,is,1,i,0))
    enddo
endif
    
print *, "kin_f"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy(0,0,is,0,i,1)),aimag(kin_energy(0,0,is,0,i,1)), real(kin_energy(0,0,is,1,i,1)),aimag(kin_energy(0,0,is,1,i,1))
    enddo



    write(j, " ")
enddo

write(j,"(a10)") "**********"

enddo


else if (l_orbital==2) then

    write(*,"(a30,i5)") "Number of n.n.:", n_nn
    write(*,"(a60)") "Nearest neighbours (c wrt to f, which is at the origin):"
    write(*,"(3a6, 5a10)") "i", "j", "k", "theta", "phi", "x", "y", "z"
    write(*,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
    write(*,"(a15)") ""

do i=1,n_nn
    write(*,"(3i6, 5f10.5)")  nn_ijk(i,1),nn_ijk(i,2),nn_ijk(i,3),nn_theta_phi(i,1),nn_theta_phi(i,2),nn_xyz(i,1),nn_xyz(i,2),nn_xyz(i,3)



    write(*,"(a35)") "Kin_energy d"
    do o1=0,1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  tc_*real(kin_energy(o1,0,is,0,i,0)),tc_*aimag(kin_energy(o1,0,is,0,i,0)), tc_*real(kin_energy(o1,0,is,1,i,0)),tc_*aimag(kin_energy(o1,0,is,1,i,0)),&
    						  tc_*real(kin_energy(o1,1,is,0,i,0)),tc_*aimag(kin_energy(o1,1,is,0,i,0)), tc_*real(kin_energy(o1,1,is,1,i,0)),tc_*aimag(kin_energy(o1,1,is,1,i,0))
    enddo
    enddo    

!print *, real(kin_energy(0,0,0,0,i,1))


    write(*,"(a35)") "Kin_energy f"
if(l_orbital_f==2) then
    do o1=0,1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  tf_*real(kin_energy(o1,0,is,0,i,1)),tf_*aimag(kin_energy(o1,0,is,0,i,1)), tf_*real(kin_energy(o1,0,is,1,i,1)),tf_*aimag(kin_energy(o1,0,is,1,i,1)),&
    						  tf_*real(kin_energy(o1,1,is,0,i,1)),tf_*aimag(kin_energy(o1,1,is,0,i,1)), tf_*real(kin_energy(o1,1,is,1,i,1)),tf_*aimag(kin_energy(o1,1,is,1,i,1))
    enddo
    enddo    
else
    do o1=0,2
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  tf_*real(kin_energy(o1,0,is,0,i,1)),tf_*aimag(kin_energy(o1,0,is,0,i,1)), tf_*real(kin_energy(o1,0,is,1,i,1)),tf_*aimag(kin_energy(o1,0,is,1,i,1)),&
    						  tf_*real(kin_energy(o1,1,is,0,i,1)),tf_*aimag(kin_energy(o1,1,is,0,i,1)), tf_*real(kin_energy(o1,1,is,1,i,1)),tf_*aimag(kin_energy(o1,1,is,1,i,1)),&
    						  tf_*real(kin_energy(o1,2,is,0,i,1)),tf_*aimag(kin_energy(o1,2,is,0,i,1)), tf_*real(kin_energy(o1,2,is,1,i,1)),tf_*aimag(kin_energy(o1,2,is,1,i,1))
    enddo
    enddo    

endif

    write(*,"(a35)") "Phi_sigma_alpha"
if(l_orbital_f==2) then
    do o1=0,1
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_sigma_alpha(o1,0,is,0,i)),Vcf_*aimag(phi_sigma_alpha(o1,0,is,0,i)), Vcf_*real(phi_sigma_alpha(o1,0,is,1,i)),Vcf_*aimag(phi_sigma_alpha(o1,0,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha(o1,1,is,0,i)),Vcf_*aimag(phi_sigma_alpha(o1,1,is,0,i)), Vcf_*real(phi_sigma_alpha(o1,1,is,1,i)),Vcf_*aimag(phi_sigma_alpha(o1,1,is,1,i))
    enddo
    enddo    
else    
    do o1=0,1
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_sigma_alpha(o1,0,is,0,i)),Vcf_*aimag(phi_sigma_alpha(o1,0,is,0,i)), Vcf_*real(phi_sigma_alpha(o1,0,is,1,i)),Vcf_*aimag(phi_sigma_alpha(o1,0,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha(o1,1,is,0,i)),Vcf_*aimag(phi_sigma_alpha(o1,1,is,0,i)), Vcf_*real(phi_sigma_alpha(o1,1,is,1,i)),Vcf_*aimag(phi_sigma_alpha(o1,1,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha(o1,2,is,0,i)),Vcf_*aimag(phi_sigma_alpha(o1,2,is,0,i)), Vcf_*real(phi_sigma_alpha(o1,2,is,1,i)),Vcf_*aimag(phi_sigma_alpha(o1,2,is,1,i))

    enddo
    enddo    
endif

    write(*,"(a35)") "Phi_alpha_sigma"

    do o1=0,l_orbital_f-1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_alpha_sigma(o1,0,is,0,i)),Vcf_*aimag(phi_alpha_sigma(o1,0,is,0,i)), Vcf_*real(phi_alpha_sigma(o1,0,is,1,i)),Vcf_*aimag(phi_alpha_sigma(o1,0,is,1,i)),&
    						  Vcf_*real(phi_alpha_sigma(o1,1,is,0,i)),Vcf_*aimag(phi_alpha_sigma(o1,1,is,0,i)), Vcf_*real(phi_alpha_sigma(o1,1,is,1,i)),Vcf_*aimag(phi_alpha_sigma(o1,1,is,1,i))
    enddo
    enddo    

    
    						 
    write(*, " ")
    enddo




endif


end subroutine write_nn_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_nn2_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)
!real(kind=idp)		:: nn_theta_phi(N_nn,2),nn_xyz(N_nn,3)
!integer, intent(in)	:: nn_ijk(N_nn,3)
!complex(kind=idpc)	::phi_sigma_alpha(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn), phi_alpha_sigma(0:l_orbital-1,0:l_orbital-1,0:1,0:1, n_nn)
!complex(kind=idpc)	::kin_energy(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn,0:1)
integer	:: i,j


if (l_orbital==1 .or. l_orbital==0) then
do j=6,24,18

    write(j,"(a30,i5)") "Number of 2nd n.n.:", n_nn2
    write(j,"(a60)") "Nearest neighbours (c wrt to f, which is at the origin):"
    write(j,"(3a6, 5a10)") "i", "j", "k", "theta", "phi", "x", "y", "z"
    write(j,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
    write(j,"(a15)") ""

do, i=1,n_nn2
    write(j,"(3i6, 5f10.5)")  nn2_ijk(i,1),nn2_ijk(i,2),nn2_ijk(i,3),nn2_theta_phi(i,1),nn2_theta_phi(i,2),nn2_xyz(i,1),nn2_xyz(i,2),nn2_xyz(i,3)

if (l_orbital>0) then
print *, "phi_sa"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_sigma_alpha2(0,0,is,0,i)),aimag(phi_sigma_alpha2(0,0,is,0,i)), real(phi_sigma_alpha2(0,0,is,1,i)),aimag(phi_sigma_alpha2(0,0,is,1,i))
    enddo

print *, "phi_as"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_alpha_sigma2(0,0,is,0,i)),aimag(phi_alpha_sigma2(0,0,is,0,i)), real(phi_alpha_sigma2(0,0,is,1,i)),aimag(phi_alpha_sigma2(0,0,is,1,i))
    enddo
print *, "kin_s"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy2(0,0,is,0,i,0)),aimag(kin_energy2(0,0,is,0,i,0)), real(kin_energy2(0,0,is,1,i,0)),aimag(kin_energy2(0,0,is,1,i,0))
    enddo
endif
    
print *, "kin_f"
    do is=0,1
    write(j,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy2(0,0,is,0,i,1)),aimag(kin_energy2(0,0,is,0,i,1)), real(kin_energy2(0,0,is,1,i,1)),aimag(kin_energy2(0,0,is,1,i,1))
    enddo




    write(j, " ")
enddo

write(j,"(a10)") "**********"

enddo


else if (l_orbital==2) then
    write(*,"(a30,i5)") "Number of 2nd n.n.:", n_nn2
    write(*,"(a60)") "Nearest neighbours (c wrt to f, which is at the origin):"
    write(*,"(3a6, 5a10)") "i", "j", "k", "theta", "phi", "x", "y", "z"
    write(*,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
    write(*,"(a15)") ""

do i=1,n_nn2
    write(*,"(4i6, 5f10.5)")  i, nn2_ijk(i,1),nn2_ijk(i,2),nn2_ijk(i,3),nn2_theta_phi(i,1),nn2_theta_phi(i,2),nn2_xyz(i,1),nn2_xyz(i,2),nn2_xyz(i,3)



    write(*,"(a35)") "Kin_energy 2 d"
    do o1=0,1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  tc_*real(kin_energy2(o1,0,is,0,i,0)),tc_*aimag(kin_energy2(o1,0,is,0,i,0)), tc_*real(kin_energy2(o1,0,is,1,i,0)),tc_*aimag(kin_energy2(o1,0,is,1,i,0)),&
    						  tc_*real(kin_energy2(o1,1,is,0,i,0)),tc_*aimag(kin_energy2(o1,1,is,0,i,0)), tc_*real(kin_energy2(o1,1,is,1,i,0)),tc_*aimag(kin_energy2(o1,1,is,1,i,0))
    enddo
    enddo    

!print *, real(kin_energy(0,0,0,0,i,1))

    write(*,"(a35)") "Kin_energy 2 f"
if (l_orbital_f==2) then
    do o1=0,1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  tf_*real(kin_energy2(o1,0,is,0,i,1)),tf_*aimag(kin_energy2(o1,0,is,0,i,1)), tf_*real(kin_energy2(o1,0,is,1,i,1)),tf_*aimag(kin_energy2(o1,0,is,1,i,1)),&
    						  tf_*real(kin_energy2(o1,1,is,0,i,1)),tf_*aimag(kin_energy2(o1,1,is,0,i,1)), tf_*real(kin_energy2(o1,1,is,1,i,1)),tf_*aimag(kin_energy2(o1,1,is,1,i,1))
    enddo
    enddo    
else if (l_orbital_f==3) then
    do o1=0,2
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  tf_*real(kin_energy2(o1,0,is,0,i,1)),tf_*aimag(kin_energy2(o1,0,is,0,i,1)), tf_*real(kin_energy2(o1,0,is,1,i,1)),tf_*aimag(kin_energy2(o1,0,is,1,i,1)),&
    						  tf_*real(kin_energy2(o1,1,is,0,i,1)),tf_*aimag(kin_energy2(o1,1,is,0,i,1)), tf_*real(kin_energy2(o1,1,is,1,i,1)),tf_*aimag(kin_energy2(o1,1,is,1,i,1)),&
    						  tf_*real(kin_energy2(o1,2,is,0,i,1)),tf_*aimag(kin_energy2(o1,2,is,0,i,1)), tf_*real(kin_energy2(o1,2,is,1,i,1)),tf_*aimag(kin_energy2(o1,2,is,1,i,1))
    enddo
    enddo    



endif





    write(*,"(a35)") "Phi_sigma_alpha"
if(l_orbital_f==2) then
    do o1=0,1
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_sigma_alpha2(o1,0,is,0,i)),Vcf_*aimag(phi_sigma_alpha2(o1,0,is,0,i)), Vcf_*real(phi_sigma_alpha2(o1,0,is,1,i)),Vcf_*aimag(phi_sigma_alpha2(o1,0,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha2(o1,1,is,0,i)),Vcf_*aimag(phi_sigma_alpha2(o1,1,is,0,i)), Vcf_*real(phi_sigma_alpha2(o1,1,is,1,i)),Vcf_*aimag(phi_sigma_alpha2(o1,1,is,1,i))
    enddo
    enddo    
else    
    do o1=0,1
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_sigma_alpha2(o1,0,is,0,i)),Vcf_*aimag(phi_sigma_alpha2(o1,0,is,0,i)), Vcf_*real(phi_sigma_alpha2(o1,0,is,1,i)),Vcf_*aimag(phi_sigma_alpha2(o1,0,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha2(o1,1,is,0,i)),Vcf_*aimag(phi_sigma_alpha2(o1,1,is,0,i)), Vcf_*real(phi_sigma_alpha2(o1,1,is,1,i)),Vcf_*aimag(phi_sigma_alpha2(o1,1,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha2(o1,2,is,0,i)),Vcf_*aimag(phi_sigma_alpha2(o1,2,is,0,i)), Vcf_*real(phi_sigma_alpha2(o1,2,is,1,i)),Vcf_*aimag(phi_sigma_alpha2(o1,2,is,1,i))

    enddo
    enddo    
endif

    write(*,"(a35)") "Phi_alpha_sigma"

    do o1=0,l_orbital_f-1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_alpha_sigma2(o1,0,is,0,i)),Vcf_*aimag(phi_alpha_sigma2(o1,0,is,0,i)), Vcf_*real(phi_alpha_sigma2(o1,0,is,1,i)),Vcf_*aimag(phi_alpha_sigma2(o1,0,is,1,i)),&
    						  Vcf_*real(phi_alpha_sigma2(o1,1,is,0,i)),Vcf_*aimag(phi_alpha_sigma2(o1,1,is,0,i)), Vcf_*real(phi_alpha_sigma2(o1,1,is,1,i)),Vcf_*aimag(phi_alpha_sigma2(o1,1,is,1,i))
    enddo
    enddo    

    
    						 
    write(*, " ")



enddo




endif


end subroutine write_nn2_phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_nn3_phi!(nn_ijk,nn_theta_phi,nn_xyz,phi_sigma_alpha, phi_alpha_sigma, kin_energy)
!real(kind=idp)		:: nn_theta_phi(N_nn,2),nn_xyz(N_nn,3)
!integer, intent(in)	:: nn_ijk(N_nn,3)
!complex(kind=idpc)	::phi_sigma_alpha(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn), phi_alpha_sigma(0:l_orbital-1,0:l_orbital-1,0:1,0:1, n_nn)
!complex(kind=idpc)	::kin_energy(0:l_orbital-1,0:l_orbital-1, 0:1,0:1, n_nn,0:1)
integer	:: i,j

    write(*,"(a30,i5)") "Number of 3rd n.n.:", n_nn3

do i=1, n_nn3
    write(*,"(4i6, 5f10.5)")  i, nn3_ijk(i,1),nn3_ijk(i,2),nn3_ijk(i,3),nn3_theta_phi(i,1),nn3_theta_phi(i,2),nn3_xyz(i,1),nn3_xyz(i,2),nn3_xyz(i,3)



if (l_orbital==1 .or. l_orbital==0) then
if (l_orbital>0) then
print *, "phi_sa"
    do is=0,1
    write(*,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_sigma_alpha3(0,0,is,0,i)),aimag(phi_sigma_alpha3(0,0,is,0,i)), real(phi_sigma_alpha3(0,0,is,1,i)),aimag(phi_sigma_alpha3(0,0,is,1,i))
    enddo

print *, "phi_as"
    do is=0,1
    write(*,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(phi_alpha_sigma3(0,0,is,0,i)),aimag(phi_alpha_sigma3(0,0,is,0,i)), real(phi_alpha_sigma3(0,0,is,1,i)),aimag(phi_alpha_sigma3(0,0,is,1,i))
    enddo
print *, "kin_s"
    do is=0,1
    write(*,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy3(0,0,is,0,i,0)),aimag(kin_energy3(0,0,is,0,i,0)), real(kin_energy3(0,0,is,1,i,0)),aimag(kin_energy3(0,0,is,1,i,0))
    enddo
endif
print *, "kin_f"
    do is=0,1
    write(*,"(2 ('(',f8.4,',',f8.4,')','  ') )")  real(kin_energy3(0,0,is,0,i,1)),aimag(kin_energy3(0,0,is,0,i,1)), real(kin_energy3(0,0,is,1,i,1)),aimag(kin_energy3(0,0,is,1,i,1))
    enddo






endif



if (l_orbital==2) then
    write(*,"(a35)") "Kin_energy 3 d"


    do o1=0,1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  tc_*real(kin_energy3(o1,0,is,0,i,0)),tc_*aimag(kin_energy3(o1,0,is,0,i,0)), tc_*real(kin_energy3(o1,0,is,1,i,0)),tc_*aimag(kin_energy3(o1,0,is,1,i,0)),&
    						  tc_*real(kin_energy3(o1,1,is,0,i,0)),tc_*aimag(kin_energy3(o1,1,is,0,i,0)), tc_*real(kin_energy3(o1,1,is,1,i,0)),tc_*aimag(kin_energy3(o1,1,is,1,i,0))
    enddo
    enddo    
endif



if (l_orbital_f==3) then
    write(*,"(a35)") "Kin_energy 3 f"

    do o1=0,2
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  tf_*real(kin_energy3(o1,0,is,0,i,1)),tf_*aimag(kin_energy3(o1,0,is,0,i,1)), tf_*real(kin_energy3(o1,0,is,1,i,1)),tf_*aimag(kin_energy3(o1,0,is,1,i,1)),&
    						  tf_*real(kin_energy3(o1,1,is,0,i,1)),tf_*aimag(kin_energy3(o1,1,is,0,i,1)), tf_*real(kin_energy3(o1,1,is,1,i,1)),tf_*aimag(kin_energy3(o1,1,is,1,i,1)),&
    						  tf_*real(kin_energy3(o1,2,is,0,i,1)),tf_*aimag(kin_energy3(o1,2,is,0,i,1)), tf_*real(kin_energy3(o1,2,is,1,i,1)),tf_*aimag(kin_energy3(o1,2,is,1,i,1))
    enddo
    enddo    
endif



    write(*,"(a35)") "Phi_sigma_alpha3"
if(l_orbital_f==3) then

    do o1=0,1
    do is=0,1
    write(*,"(6 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_sigma_alpha3(o1,0,is,0,i)),Vcf_*aimag(phi_sigma_alpha3(o1,0,is,0,i)), Vcf_*real(phi_sigma_alpha3(o1,0,is,1,i)),Vcf_*aimag(phi_sigma_alpha3(o1,0,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha3(o1,1,is,0,i)),Vcf_*aimag(phi_sigma_alpha3(o1,1,is,0,i)), Vcf_*real(phi_sigma_alpha3(o1,1,is,1,i)),Vcf_*aimag(phi_sigma_alpha3(o1,1,is,1,i)),&
    						  Vcf_*real(phi_sigma_alpha3(o1,2,is,0,i)),Vcf_*aimag(phi_sigma_alpha3(o1,2,is,0,i)), Vcf_*real(phi_sigma_alpha3(o1,2,is,1,i)),Vcf_*aimag(phi_sigma_alpha3(o1,2,is,1,i))

    enddo
    enddo    

    write(*,"(a35)") "Phi_alpha_sigma3"

    do o1=0,l_orbital_f-1
    do is=0,1
    write(*,"(4 ('(',f8.4,',',f8.4,')','  ') )")  Vcf_*real(phi_alpha_sigma3(o1,0,is,0,i)),Vcf_*aimag(phi_alpha_sigma3(o1,0,is,0,i)), Vcf_*real(phi_alpha_sigma3(o1,0,is,1,i)),Vcf_*aimag(phi_alpha_sigma3(o1,0,is,1,i)),&
    						  Vcf_*real(phi_alpha_sigma3(o1,1,is,0,i)),Vcf_*aimag(phi_alpha_sigma3(o1,1,is,0,i)), Vcf_*real(phi_alpha_sigma3(o1,1,is,1,i)),Vcf_*aimag(phi_alpha_sigma3(o1,1,is,1,i))
    enddo
    enddo    

endif

    write(*, " ")

enddo

end subroutine write_nn3_phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_states

 write(*,"(a20)")  "List of states:"
 write(*,"(10a6)")  "ind", "x", "y", "z", "c/f", "s"
 do indp=0,total_size-1
 call convert_index(indp, ix, iy, iz, icf, is)
 write(*,"(10i6)")  indp, ix, iy, iz, icf, is
 enddo
 
 write(*,"(a20)")  "List of sites:"
 write(*,"(10a6)")  "ind", "x", "y", "z"
 do ind=0,n_sites-1
 call convert_index_site(ind, ix, iy, iz)
 write(*,"(10i6)")  ind, ix, iy, iz
 enddo

write(*,"(a10)") "**********"

end subroutine write_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_cf_matrix(alpha_jzi_mat)
real(kind=idp)		:: alpha_jzi_mat(0:5,0:5)

alpha_jzi_mat=0

if (orb_type=="f5")  then

!2nd index: 0:-5-2   1:-3-2   2:-1-2   3:1-2   4:3-2   5:5-2
alpha_jzi_mat(0,0)=+sqrt(5.0d0/6.0d0)	!gamma8(1)-
alpha_jzi_mat(0,4)=+sqrt(1.0d0/6.0d0)
alpha_jzi_mat(1,1)=+sqrt(1.0d0/6.0d0)	!gamma8(1)+
alpha_jzi_mat(1,5)=+sqrt(5.0d0/6.0d0)
alpha_jzi_mat(2,2)=1.0d0		!gamma8(2)-
alpha_jzi_mat(3,3)=1.0d0		!gamma8(2)+	
alpha_jzi_mat(4,0)=+sqrt(1.0d0/6.0d0)	!gamma7-
alpha_jzi_mat(4,4)=-sqrt(5.0d0/6.0d0)
alpha_jzi_mat(5,1)=-sqrt(5.0d0/6.0d0)	!gamma7+
alpha_jzi_mat(5,5)=+sqrt(1.0d0/6.0d0)

    write(*,"(a30)") ( "Crystal field matrix:" )
    write(*,"(a5,10a10)") "", "-5-2","-3-2", "-1-2", "1-2", "3-2", "5-2"

    write(*,"(a5,10(f10.5))") "g1-", ( alpha_jzi_mat(0,j) , j=0,5 )
    write(*,"(a5,10(f10.5))") "g1+", ( alpha_jzi_mat(1,j) , j=0,5 )
    write(*,"(a5,10(f10.5))") "g2-", ( alpha_jzi_mat(2,j) , j=0,5 )
    write(*,"(a5,10(f10.5))") "g2+", ( alpha_jzi_mat(3,j) , j=0,5 )
    write(*,"(a5,10(f10.5))") "g3-", ( alpha_jzi_mat(4,j) , j=0,5 )
    write(*,"(a5,10(f10.5))") "g3+", ( alpha_jzi_mat(5,j) , j=0,5 )


elseif (orb_type=="p3") then !J=3/2
alpha_jzi_mat(0,1)=1.0d0	!-3/2
alpha_jzi_mat(1,4)=1.0d0	!+3/2
alpha_jzi_mat(2,2)=1.0d0	!-1/2
alpha_jzi_mat(3,3)=1.0d0	!+1/2

elseif (orb_type=="p1") then !J=1/2
alpha_jzi_mat(0,2)=1.0d0	!-1/2
alpha_jzi_mat(1,3)=1.0d0	!+1/2

endif

UtU=0

do i=0,5
 do j=0,5
  do k=0,5
  UtU(i,j)=UtU(i,j)+alpha_jzi_mat(k,i)*alpha_jzi_mat(k,j)
  enddo
 enddo
enddo

   write(*,"(a30)") ( "UtU:" )
 do i=0,5
    write(*,"(10(f10.5))") ( UtU(i,j) , j=0,5 )
enddo

end subroutine set_cf_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_doublet!(alpha_jzi, alpha_msigma)
!real(kind=idp)		::alpha_jzi(0:l_orbital-1,0:1,0:5), alpha_msigma(0:l_orbital-1,0:1,0:1,-3:3), norm(0:l_orbital-1,0:1)
real(kind=idp)		::norm(0:l_orbital_f-1,0:1)

alpha_jzi=0
!standard: Gamma_8(2)
alpha_jzi(0,0,:)=alpha_jzi_mat(2,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(3,:)

!!!!!!!!!Gamma_1: Gamma_8(1)
if (phi_type=="g_1") then
alpha_jzi(0,0,:)=alpha_jzi_mat(0,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(1,:)

!!!!!!!!!Gamma_2: Gamma_8(2)
elseif (phi_type=="g_2") then
alpha_jzi(0,0,:)=alpha_jzi_mat(2,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(3,:)

!!!!!!!!!Gamma_3 // Gamma_7
elseif (phi_type=="g_3") then
alpha_jzi(0,0,:)=alpha_jzi_mat(4,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(5,:)

!!!!!!!!!Jz=+-1-2
elseif (phi_type=="1-2" .or. phi_type=="sxy") then
alpha_jzi(0,0,2)=1	!-|-1/2>
alpha_jzi(0,1,3)=1	!+|+1/2>

!!!!!!!!!Jz=+-3-2
elseif (phi_type=="3-2") then
alpha_jzi(0,0,1)=-1
alpha_jzi(0,1,4)=1

!!!!!!!!!Jz=+-5-2
elseif (phi_type=="5-2") then
alpha_jzi(0,0,0)=-1
alpha_jzi(0,1,5)=1

!elseif (phi_type=="g7_".or. phi_type=="sm2") then
elseif (phi_type=="g7_") then
alpha_jzi(0,0,:)=alpha_jzi_mat(4,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(5,:)




!0-1:gamma8(1) 2-3 gamma8(2) 4-5:gamma7

elseif (phi_type=="smb") then
alpha_jzi(0,0,:)=alpha_jzi_mat(0,:)	!f1- Gamma_8(1)
alpha_jzi(0,1,:)=alpha_jzi_mat(1,:)	!f1+

alpha_jzi(1,0,:)=alpha_jzi_mat(2,:)	!f2-=|-1/2> Gamma_8(2)
alpha_jzi(1,1,:)=alpha_jzi_mat(3,:)	!f2+=|+1/2>

elseif (phi_type=="sm7") then
alpha_jzi(0,0,:)=alpha_jzi_mat(0,:)	!f1- Gamma_81-
alpha_jzi(0,1,:)=alpha_jzi_mat(1,:)	!f1+ Gamma_81+
alpha_jzi(1,0,:)=alpha_jzi_mat(2,:)	!f2- Gamma_82-
alpha_jzi(1,1,:)=alpha_jzi_mat(3,:)	!f2+ Gamma_82+
alpha_jzi(2,0,:)=alpha_jzi_mat(4,:)	!f3- Gamma_7-
alpha_jzi(2,1,:)=alpha_jzi_mat(5,:)	!f3+ Gamma_7+


!p states
elseif (phi_type=="p33") then	!3/2 3/2
alpha_jzi(0,0,:)=alpha_jzi_mat(0,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(1,:)

elseif (phi_type=="p31"  .or. phi_type=="pd4") then	!3/2 1/2
alpha_jzi(0,0,:)=alpha_jzi_mat(2,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(3,:)

elseif (phi_type=="p11"  .or. phi_type=="pd2") then	!1/2 1/2
alpha_jzi(0,0,:)=alpha_jzi_mat(0,:)
alpha_jzi(0,1,:)=alpha_jzi_mat(1,:)



endif


write(*,"(a10)") "**********"
print *, "Alfa doublet over J_z (from input):"
write(*,"(100a6)") "","-5-2", "-3-2", "-1-2", "1-2", "3-2", "5-2"

if (l_orbital_f==1) then
    write(*,"(a6,100f6.2)"),"-" ,( alpha_jzi(0,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"+", ( alpha_jzi(0,1,j), j=0,5 )
else if (l_orbital_f==2) then
    write(*,"(a6,100f6.2)"),"1-" ,( alpha_jzi(0,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"1+", ( alpha_jzi(0,1,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"2-" ,( alpha_jzi(1,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"2+", ( alpha_jzi(1,1,j), j=0,5 )
else if (l_orbital_f==3) then
    write(*,"(a6,100f6.2)"),"1-" ,( alpha_jzi(0,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"1+", ( alpha_jzi(0,1,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"2-" ,( alpha_jzi(1,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"2+", ( alpha_jzi(1,1,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"3-" ,( alpha_jzi(2,0,j), j=0,5 )
    write(*,"(a6,100f6.2)"),"3+", ( alpha_jzi(2,1,j), j=0,5 )
endif

alpha_msigma=0	!alpha state over sigma and m
do i=0,1	!alpha
do o1=0,l_orbital_f-1
 do j=0,1	!sigma
  do k=-3,3	!m
   do l=0,5	!J_z
   alpha_msigma(o1,i,j,k)=alpha_msigma(o1,i,j,k)+alpha_jzi(o1,i,l)*cg(dble(j1),0.5d0,j3,dble(k),dble(j)-0.5d0,dble(l)-2.5d0)
   enddo
  enddo
 enddo
enddo
enddo

norm=0
do l=0,1
do o1=0,l_orbital_f-1
 do j=0,1	!sigma
  do k=-3,3	!m
  norm(o1,l)=norm(o1,l)+(abs(alpha_msigma(o1,l,j,k)))**2
  enddo
 enddo
enddo
enddo

    write(*,"(a100)") "Alfa doublet (+/-) over sigma(up/down),m(-3->+3):"
    write(*,"(a6,7a10)"),"","-3","-2","-1","0","+1","+2","+3"
if (l_orbital_f==1) then
    write(*,"(a6,100f10.5)"),"- down", ( alpha_msigma(0,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-   up", ( alpha_msigma(0,0,1,k), k=-3,3 ), norm(0,0)
    write(*,"(a6,100f10.5)"),"+ down", ( alpha_msigma(0,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+   up", ( alpha_msigma(0,1,1,k), k=-3,3 ), norm(0,1)


elseif (l_orbital_f==2) then
    write(*,"(a6,100f10.5)"),"-1 down", ( alpha_msigma(0,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-1   up", ( alpha_msigma(0,0,1,k), k=-3,3 ), norm(0,0)
    write(*,"(a6,100f10.5)"),"+1 down", ( alpha_msigma(0,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+1  up", ( alpha_msigma(0,1,1,k), k=-3,3 ), norm(0,1)
    write(*,"(a6,100f10.5)"),"-2 down", ( alpha_msigma(1,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-2   up", ( alpha_msigma(1,0,1,k), k=-3,3 ), norm(1,0)
    write(*,"(a6,100f10.5)"),"+2 down", ( alpha_msigma(1,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+2   up", ( alpha_msigma(1,1,1,k), k=-3,3 ), norm(1,1)

elseif (l_orbital_f==3) then
    write(*,"(a6,100f10.5)"),"-1 down", ( alpha_msigma(0,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-1   up", ( alpha_msigma(0,0,1,k), k=-3,3 ), norm(0,0)
    write(*,"(a6,100f10.5)"),"+1 down", ( alpha_msigma(0,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+1  up", ( alpha_msigma(0,1,1,k), k=-3,3 ), norm(0,1)
    write(*,"(a6,100f10.5)"),"-2 down", ( alpha_msigma(1,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-2   up", ( alpha_msigma(1,0,1,k), k=-3,3 ), norm(1,0)
    write(*,"(a6,100f10.5)"),"+2 down", ( alpha_msigma(1,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+2   up", ( alpha_msigma(1,1,1,k), k=-3,3 ), norm(1,1)
    write(*,"(a6,100f10.5)"),"-3 down", ( alpha_msigma(2,0,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"-3   up", ( alpha_msigma(2,0,1,k), k=-3,3 ), norm(2,0)
    write(*,"(a6,100f10.5)"),"+3 down", ( alpha_msigma(2,1,0,k), k=-3,3 )
    write(*,"(a6,100f10.5)"),"+3   up", ( alpha_msigma(2,1,1,k), k=-3,3 ), norm(2,1)


endif

    write(*,"(a10)") "**********"



    write(24,"(a100)") "Alfa doublet (+/-) over sigma(up/down),m(-3->+3):"
    write(24,"(a6,7a10)"),"","-3","-2","-1","0","+1","+2","+3"
    write(24,"(a6,100f10.5)"),"- down", ( alpha_msigma(0,0,0,k), k=-3,3 )
    write(24,"(a6,100f10.5)"),"-   up", ( alpha_msigma(0,0,1,k), k=-3,3 ), norm(0,0)
    write(24,"(a6,100f10.5)"),"+ down", ( alpha_msigma(0,1,0,k), k=-3,3 )
    write(24,"(a6,100f10.5)"),"+   up", ( alpha_msigma(0,1,1,k), k=-3,3 ), norm(0,1)
    write(24,"(a10)") "**********"

!print *, alpha_msigma


end subroutine set_doublet


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rotation in k-space: must have nx=nkx, ny=nky, nz=nkz, n_sites=n_k
subroutine set_U(U, Um1)
complex(kind=idpc), intent(out)	::U(0:total_size-1,0:total_size-1),Um1(0:total_size-1,0:total_size-1)
complex(kind=idpc)	::UstarU(0:total_size-1,0:total_size-1)
integer		::ind,ix,iy,iz,ik, i,j,k


if (n_k .ne. n_sites) then
print *, "n_k .ne. n_sites"
stop
endif

U=0

!f_k_c=f_k_c+ exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs(ik,1)*lkx+nn_ijk(i,2)*k_vecs(ik,2)*lky+nn_ijk(i,3)*k_vecs(ik,3)*lkz))


do ind=0,n_sites-1
  do ik=0,n_k-1

   call convert_index_site(ind,ix,iy,iz)
 
! U(ind,ik)= exp(2*(0,1)*pi*(ix*k_vecs(ik,1)+iy*k_vecs(ik,2)+iz*k_vecs(ik,3)))/sqrt(real(n_sites))				!cup
 U(ind+n_sites,ik+n_k)=U(ind,ik)		!fup
 U(ind+2*n_sites,ik+2*n_k)=U(ind,ik)		!cdown
 U(ind+3*n_sites,ik+3*n_k)=U(ind,ik)		!fdown

  enddo
enddo


do i=0,total_size-1
  do j=0,total_size-1

   Um1(i,j)=conjg(U(j,i))
 
  enddo
enddo


UstarU=0

do i=0,total_size-1
  do j=0,total_size-1
     do k=0,total_size-1

   UstarU (i,j)=UstarU(i,j)+Um1(i,k)*U(k,j)
   
     enddo
  enddo
enddo
	
!if (n_sites < 5) call write_ham(UstarU,"UstarU")


do i=0,total_size-1
  do j=0,total_size-1

   if (i==j .and. abs(Ustaru(i,j)-1)>1e-6) then
   print *, "error in U",i, j , Ustaru(i,j)
 !  stop
   endif

   if (i.ne.j .and. abs(Ustaru(i,j))>1e-6) then
   print *, "error in U", i, j , Ustaru(i,j)
!   stop
   endif

  enddo
enddo



	
end subroutine set_U



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BASIS_ROTATIONcmplx (U2_matrix, A_matrix, U1_matrix, result_matrix)

IMPLICIT NONE

COMPLEX(KIND=idpc), INTENT(IN)     :: A_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(IN)     :: U1_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(IN)     :: U2_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(INOUT)     :: result_matrix(:,:)
COMPLEX(KIND=idpc), ALLOCATABLE :: tmp_matrix(:,:)

!result_matrix = MATMUL( TRANSPOSE(real(U2_matrix)-(0,1)*aimag(U2_matrix)), MATMUL(A_matrix, U1_matrix) )
!result_matrix = MATMUL( TRANSPOSE(U2_matrix), MATMUL(A_matrix, U1_matrix) )



! A * U1
ALLOCATE( tmp_matrix( SIZE(A_matrix, 1), SIZE(U1_matrix, 2) ) )
tmp_matrix = 0.0D0

CALL ZGEMM( 'N', 'N', SIZE(A_matrix, 1), SIZE(U1_matrix, 2), &
        SIZE(A_matrix, 2), 1.0D0, A_matrix, SIZE(A_matrix, 1), &
        U1_matrix, SIZE(U1_matrix, 1), 0.0D0, &
        tmp_matrix, SIZE(tmp_matrix, 1) )

! U2^+ * tmp
result_matrix=0

CALL ZGEMM( 'N', 'N', SIZE(U2_matrix, 2), SIZE(tmp_matrix, 2), &
        SIZE(U2_matrix, 1), 1.0D0, U2_matrix, SIZE(U2_matrix, 1), &
       tmp_matrix, SIZE(tmp_matrix, 1), 0.0D0, &
        result_matrix, SIZE(result_matrix, 1) )!

DEALLOCATE( tmp_matrix )




END SUBROUTINE BASIS_ROTATIONcmplx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!transpose conjugate (U2)*A*U1
!A(N,N)
!Result(M,M)
SUBROUTINE BASIS_ROTATIONcmplx_c (U2_matrix, A_matrix, U1_matrix, result_matrix,n,m)

IMPLICIT NONE

integer				   ::M,N
COMPLEX(KIND=idpc), INTENT(IN)     :: A_matrix(N,N)
COMPLEX(KIND=idpc), INTENT(IN)     :: U1_matrix(N,M)
COMPLEX(KIND=idpc), INTENT(IN)     :: U2_matrix(N,M)
COMPLEX(KIND=idpc), INTENT(INOUT)     :: result_matrix(M,M)
COMPLEX(KIND=idpc), ALLOCATABLE :: tmp_matrix(:,:)


! A * U1
ALLOCATE( tmp_matrix( N,M ))
tmp_matrix = 0.0D0

CALL ZGEMM( 'N', 'N', N, M, N, 1.0D0, A_matrix, N, U1_matrix, N, 0.0D0, tmp_matrix, N )

! U2^+ * tmp
result_matrix=0

CALL ZGEMM( 'C', 'N', M, M, N, 1.0D0, U2_matrix, N,tmp_matrix, N, 0.0D0, result_matrix, M )!

DEALLOCATE( tmp_matrix )

result_matrix = MATMUL( TRANSPOSE(real(U2_matrix)-(0,1)*aimag(U2_matrix)), MATMUL(A_matrix, U1_matrix) )


END SUBROUTINE BASIS_ROTATIONcmplx_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BASIS_ROTATIONcmplxhc (U2_matrix, A_matrix, U1_matrix, result_matrix)

IMPLICIT NONE

COMPLEX(KIND=idpc), INTENT(IN)     :: A_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(IN)     :: U1_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(IN)     :: U2_matrix(:,:)
COMPLEX(KIND=idpc), INTENT(INOUT)     :: result_matrix(:,:)
COMPLEX(KIND=idpc), ALLOCATABLE :: tmp_matrix(:,:)

result_matrix = MATMUL( U2_matrix, MATMUL(A_matrix, TRANSPOSE(real(U1_matrix)-(0,1)*aimag(U1_matrix))) )
!result_matrix = MATMUL( TRANSPOSE(real(U2_matrix)-(0,1)*aimag(U2_matrix)), MATMUL(A_matrix, U1_matrix) )
!result_matrix = MATMUL( TRANSPOSE(U2_matrix), MATMUL(A_matrix, U1_matrix) )



! A * U1
!ALLOCATE( tmp_matrix( SIZE(A_matrix, 1), SIZE(U1_matrix, 2) ) )
!tmp_matrix = 0.0D0

!CALL ZGEMM( 'N', 'C', SIZE(A_matrix, 1), SIZE(U1_matrix, 2), &
!        SIZE(A_matrix, 2), 1.0D0, A_matrix, SIZE(A_matrix, 1), &
!        U1_matrix, SIZE(U1_matrix, 1), 0.0D0, &
!        tmp_matrix, SIZE(tmp_matrix, 1) )

! U2^+ * tmp
!result_matrix=0

!CALL ZGEMM( 'N', 'N', SIZE(U2_matrix, 2), SIZE(tmp_matrix, 2), &
!        SIZE(U2_matrix, 1), 1.0D0, U2_matrix, SIZE(U2_matrix, 1), &
!       tmp_matrix, SIZE(tmp_matrix, 1), 0.0D0, &
!        result_matrix, SIZE(result_matrix, 1) )!

!DEALLOCATE( tmp_matrix )

END SUBROUTINE BASIS_ROTATIONcmplxhc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,b,lambda,mu,ldiag,nx,ny,nz,ec_site,ef_site, k_vecs_r, n_k)
integer, intent(in)	::ix,iy,iz,icf,is,ik,o1,nx,ny,n_k,nz
integer, intent(out)	::index_psi(n_max_link,8)
complex(kind=idpc), intent(out)	::h_psi(n_max_link)
logical,intent(in)	::ldiag
integer	::deltanx,deltany,deltanz

!real(kind=idp), intent(in)	:: b(0:n_sites-1), lambda(0:n_sites-1),mu
real(kind=idp), intent(in)	:: b(0:nx*ny*nz-1), lambda(0:nx*ny*nz-1),mu
real (kind=idp),intent(in)	:: ec_site(0:nx*ny*nz-1,0:max(1,l_orbital)-1), ef_site(0:nx*ny*nz-1,0:l_orbital_f-1), k_vecs_r(0:n_k-1,4)



ind=index_xyz(ix,iy,iz,nx,ny)
   

index_psi=0
h_psi=0

if (.not. lreadham) then
!!!!!!!!!on site: atm only energy 1->dim_hk
j=1

ind2=ind

!do icf2=0,1
!do is2=0,1
!do o2=0,l_orb(icf2)-1


 index_psi(j,1)=ix
 index_psi(j,2)=iy
 index_psi(j,3)=iz
 index_psi(j,4)=icf
 index_psi(j,5)=is
 index_psi(j,8)=o1

if(ldiag) then
if (icf==0) H_psi(j)=ec_site(ind,o1)-mu
if (icf==1) H_psi(j)=ef_site(ind,o1)-lambda(ind)-mu	
endif



!onsite hopping
j=j+1
 index_psi(j,1)=ix
 index_psi(j,2)=iy
 index_psi(j,3)=iz
 index_psi(j,4)=1-icf
 index_psi(j,5)=is
 index_psi(j,8)=o1


if (phi_type=="ons") then
H_psi(j)=Vcf_*b(ind)
endif


!enddo
!enddo
!enddo


!!!!!!!!!!!
!first neighbors
do i=1,n_nn
 ix2=mod(ix-(nn_ss(1,1)*nn_ijk(i,1)+nn_ss(1,2)*nn_ijk(i,2)+nn_ss(1,3)*nn_ijk(i,3))+9*nx,nx)		!r=r'+l
 iy2=mod(iy-(nn_ss(2,1)*nn_ijk(i,1)+nn_ss(2,2)*nn_ijk(i,2)+nn_ss(2,3)*nn_ijk(i,3))+9*ny,ny)
 iz2=mod(iz-(nn_ss(3,1)*nn_ijk(i,1)+nn_ss(3,2)*nn_ijk(i,2)+nn_ss(3,3)*nn_ijk(i,3))+9*nz,nz)  
 deltanx=-(ix-(nn_ss(1,1)*nn_ijk(i,1)+nn_ss(1,2)*nn_ijk(i,2)+nn_ss(1,3)*nn_ijk(i,3))-ix2)/nx	!how many supercells I have to cross --> to me multiplied for the exp factor with k points
 deltany=-(iy-(nn_ss(2,1)*nn_ijk(i,1)+nn_ss(2,2)*nn_ijk(i,2)+nn_ss(2,3)*nn_ijk(i,3))-iy2)/ny
 deltanz=-(iz-(nn_ss(3,1)*nn_ijk(i,1)+nn_ss(3,2)*nn_ijk(i,2)+nn_ss(3,3)*nn_ijk(i,3))-iz2)/nz

!write(*, "(60i5)"),ix,ix2,iy,iy2,iz,iz2, nn_ijk(i,1), nn_ijk(i,2), nn_ijk(i,3)
ind2=index_xyz(ix2,iy2,iz2,nx,ny)

do icf2=0,1
do is2=0,1
do o2=0,l_orb(icf2)-1
 
 !j=index_cso(icf2,is2,o2)+2+dim_hk*(i-1)
j=j+1
 
 index_psi(j,1)=ix2
 index_psi(j,2)=iy2
 index_psi(j,3)=iz2
 index_psi(j,4)=icf2
 index_psi(j,5)=is2
 index_psi(j,8)=o2

if (icf==0 .and. icf2==0) H_psi(j)=kin_energy(o1,o2,is,is2,i,0)*tc_
if (icf==1 .and. icf2==1) H_psi(j)=kin_energy(o1,o2,is,is2,i,1)*tf_*b(ind)*b(ind2)
if (icf==0 .and. icf2==1) H_psi(j)=phi_sigma_alpha(o1,o2,is,is2,i)*Vcf_*b(ind2)
if (icf==1 .and. icf2==0) H_psi(j)=phi_alpha_sigma(o1,o2,is,is2,i)*Vcf_*b(ind)


!boundaries
 if (.not. pbcx .and. deltanx.ne.0)  H_psi(j)=0
 if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanx*k_vecs_r(ik,1))
 if (.not. pbcy  .and. deltany.ne.0)  H_psi(j)=0
 if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltany*k_vecs_r(ik,2))
 if (.not. pbcz .and. deltanz.ne.0)  H_psi(j)=0
 if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanz*k_vecs_r(ik,3))
!if (ix==ix2+nn_ijk(i,1)-nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (ix==ix2+nn_ijk(i,1)+nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (iy==iy2+nn_ijk(i,2)-ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iy==iy2+nn_ijk(i,2)+ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iz==iz2+nn_ijk(i,3)-nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,3))
!endif
!if (iz==iz2+nn_ijk(i,3)+nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,3))
!endif

enddo
enddo
enddo
enddo	!i


!second neighbors
do i=1,n_nn2
 ix2=mod(ix-(nn_ss(1,1)*nn2_ijk(i,1)+nn_ss(1,2)*nn2_ijk(i,2)+nn_ss(1,3)*nn2_ijk(i,3))+9*nx,nx)		!r=r'+l
 iy2=mod(iy-(nn_ss(2,1)*nn2_ijk(i,1)+nn_ss(2,2)*nn2_ijk(i,2)+nn_ss(2,3)*nn2_ijk(i,3))+9*ny,ny)
 iz2=mod(iz-(nn_ss(3,1)*nn2_ijk(i,1)+nn_ss(3,2)*nn2_ijk(i,2)+nn_ss(3,3)*nn2_ijk(i,3))+9*nz,nz)  
ind2=index_xyz(ix2,iy2,iz2,nx,ny)
 deltanx=-(ix-(nn_ss(1,1)*nn2_ijk(i,1)+nn_ss(1,2)*nn2_ijk(i,2)+nn_ss(1,3)*nn2_ijk(i,3))-ix2)/nx	!how many supercells I have to cross --> to me multiplied for the exp factor with k points
 deltany=-(iy-(nn_ss(2,1)*nn2_ijk(i,1)+nn_ss(2,2)*nn2_ijk(i,2)+nn_ss(2,3)*nn2_ijk(i,3))-iy2)/ny
 deltanz=-(iz-(nn_ss(3,1)*nn2_ijk(i,1)+nn_ss(3,2)*nn2_ijk(i,2)+nn_ss(3,3)*nn2_ijk(i,3))-iz2)/nz

!write(*, "(60i5)"),ix,ix2,iy,iy2,iz,iz2, nn2_ijk(i,1), nn2_ijk(i,2), nn2_ijk(i,3)


do icf2=0,1
do is2=0,1
do o2=0,l_orb(icf2)-1
 
! j=index_cso(icf2,is2,o2)+2+dim_hk*(n_nn)+dim_hk*(i-1)
 j=j+1
 
 index_psi(j,1)=ix2
 index_psi(j,2)=iy2
 index_psi(j,3)=iz2
 index_psi(j,4)=icf2
 index_psi(j,5)=is2
 index_psi(j,8)=o2

if (icf==0 .and. icf2==0) H_psi(j)=kin_energy2(o1,o2,is,is2,i,0)*tc_
if (icf==1 .and. icf2==1) H_psi(j)=kin_energy2(o1,o2,is,is2,i,1)*tf_*b(ind)*b(ind2)
if (icf==0 .and. icf2==1) H_psi(j)=phi_sigma_alpha2(o1,o2,is,is2,i)*Vcf_*b(ind2)
if (icf==1 .and. icf2==0) H_psi(j)=phi_alpha_sigma2(o1,o2,is,is2,i)*Vcf_*b(ind)

!boundaries
 if (.not. pbcx .and. deltanx.ne.0)  H_psi(j)=0
 if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanx*k_vecs_r(ik,1))
 if (.not. pbcy  .and. deltany.ne.0)  H_psi(j)=0
 if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltany*k_vecs_r(ik,2))
 if (.not. pbcz .and. deltanz.ne.0)  H_psi(j)=0
 if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanz*k_vecs_r(ik,3))
!if (ix==ix2+nn2_ijk(i,1)-nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (ix==ix2+nn2_ijk(i,1)+nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (iy==iy2+nn2_ijk(i,2)-ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iy==iy2+nn2_ijk(i,2)+ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iz==iz2+nn2_ijk(i,3)-nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,3))
!endif
!if (iz==iz2+nn2_ijk(i,3)+nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,3))
!endif

enddo
enddo
enddo
enddo

!third neighbors
do i=1,n_nn3
 ix2=mod(ix-(nn_ss(1,1)*nn3_ijk(i,1)+nn_ss(1,2)*nn3_ijk(i,2)+nn_ss(1,3)*nn3_ijk(i,3))+9*nx,nx)		!r=r'+l
 iy2=mod(iy-(nn_ss(2,1)*nn3_ijk(i,1)+nn_ss(2,2)*nn3_ijk(i,2)+nn_ss(2,3)*nn3_ijk(i,3))+9*ny,ny)
 iz2=mod(iz-(nn_ss(3,1)*nn3_ijk(i,1)+nn_ss(3,2)*nn3_ijk(i,2)+nn_ss(3,3)*nn3_ijk(i,3))+9*nz,nz)  
ind2=index_xyz(ix2,iy2,iz2,nx,ny)
 deltanx=-(ix-(nn_ss(1,1)*nn3_ijk(i,1)+nn_ss(1,2)*nn3_ijk(i,2)+nn_ss(1,3)*nn3_ijk(i,3))-ix2)/nx	!how many supercells I have to cross --> to me multiplied for the exp factor with k points
 deltany=-(iy-(nn_ss(2,1)*nn3_ijk(i,1)+nn_ss(2,2)*nn3_ijk(i,2)+nn_ss(2,3)*nn3_ijk(i,3))-iy2)/ny
 deltanz=-(iz-(nn_ss(3,1)*nn3_ijk(i,1)+nn_ss(3,2)*nn3_ijk(i,2)+nn_ss(3,3)*nn3_ijk(i,3))-iz2)/nz

!write(*, "(60i5)"),ix,ix2,iy,iy2,iz,iz2, nn3_ijk(i,1), nn3_ijk(i,2), nn3_ijk(i,3)

do icf2=0,1
do is2=0,1
do o2=0,l_orb(icf2)-1

 j=j+1
! j=index_cso(icf2,is2,o2)+2+dim_hk*(n_nn+n_nn2)+dim_hk*(i-1)
 
 index_psi(j,1)=ix2
 index_psi(j,2)=iy2
 index_psi(j,3)=iz2
 index_psi(j,4)=icf2
 index_psi(j,5)=is2
 index_psi(j,8)=o2

if (icf==0 .and. icf2==0) H_psi(j)=kin_energy3(o1,o2,is,is2,i,0)*tc_
if (icf==1 .and. icf2==1) H_psi(j)=kin_energy3(o1,o2,is,is2,i,1)*tf_*b(ind)*b(ind2)
if (icf==0 .and. icf2==1) H_psi(j)=phi_sigma_alpha3(o1,o2,is,is2,i)*Vcf_*b(ind2)
if (icf==1 .and. icf2==0) H_psi(j)=phi_alpha_sigma3(o1,o2,is,is2,i)*Vcf_*b(ind)

!boundaries
 if (.not. pbcx .and. deltanx.ne.0)  H_psi(j)=0
 if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanx*k_vecs_r(ik,1))
 if (.not. pbcy  .and. deltany.ne.0)  H_psi(j)=0
 if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltany*k_vecs_r(ik,2))
 if (.not. pbcz .and. deltanz.ne.0)  H_psi(j)=0
 if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanz*k_vecs_r(ik,3))
!if (ix==ix2+nn3_ijk(i,1)-nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (ix==ix2+nn3_ijk(i,1)+nx) then
! if (.not. pbcx)  H_psi(j)=0
! if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,1))
!endif
!if (iy==iy2+nn3_ijk(i,2)-ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iy==iy2+nn3_ijk(i,2)+ny) then
! if (.not. pbcy)  H_psi(j)=0
! if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,2))
!endif
!if (iz==iz2+nn3_ijk(i,3)-nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*k_vecs_r(ik,3))
!endif
!if (iz==iz2+nn3_ijk(i,3)+nz) then
! if (.not. pbcz)  H_psi(j)=0
! if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(-2*pi*(0,1)*k_vecs_r(ik,3))
!endif


enddo
enddo
enddo
enddo

!!!!!!!!!!!!!!Reading
else !readham

ind=index_xyz(ix,iy,iz,nx,ny)
i=index_cso(icf,is,o1)

j=1

!onsite energy
if (ldiag) then
 index_psi(j,1)=ix
 index_psi(j,2)=iy
 index_psi(j,3)=iz
 index_psi(j,4)=icf
 index_psi(j,5)=is
 index_psi(j,8)=o1

 if (icf==0)  H_psi(j)=Ham_input(0,0,0,i,i)-mu+ec_site(ind,o1)			
 if (icf==1)  H_psi(j)=Ham_input(0,0,0,i,i)-mu+ef_site(ind,o1)-lambda(ind)
endif

do lx=-n_shells,n_shells
do ly=-n_shells,n_shells
do lz=-n_shells,n_shells

if (sqrt(real(lx**2+ly**2+lz**2)) <=(dist_max+0.001) .and. sqrt(real(lx**2+ly**2+lz**2)) >0.001 ) then !I exclude on site and too far shells


 ix2=mod(ix-lx+n_shells*nx,nx)	!add n_shells*nx to get >0 and get the mod correct
 iy2=mod(iy-ly+n_shells*ny,ny)
 iz2=mod(iz-lz+n_shells*nz,nz)  
 ind2=index_xyz(ix2,iy2,iz2,nx,ny)
 deltanx=-(ix-lx-ix2)/nx	!how many supercells I have to cross --> to me multiplied for the exp factor with k points
 deltany=-(iy-ly-iy2)/ny
 deltanz=-(iz-lz-iz2)/nz
 
do icf2=0,1
do is2=0,1
do o2=0,l_orb(icf2)-1
j=j+1
!write (*, "(5i5)"), j, lx,ly,lz

l=index_cso(icf2,is2,o2)

 index_psi(j,1)=ix2
 index_psi(j,2)=iy2
 index_psi(j,3)=iz2
 index_psi(j,4)=icf2
 index_psi(j,5)=is2
 index_psi(j,8)=o2


!dd
if (icf==0 .and. icf2==0) H_psi(j)=Ham_input(lx,ly,lz,i,l)
!df
if (icf==1 .and. icf2==0) H_psi(j)=Ham_input(lx,ly,lz,i,l)*b(ind)
if (icf==0 .and. icf2==1) H_psi(j)=Ham_input(lx,ly,lz,i,l)*b(ind2)
!ff
if (icf==1 .and. icf2==1) H_psi(j)=Ham_input(lx,ly,lz,i,l)*b(ind)*b(ind2)


 if (.not. pbcx .and. deltanx.ne.0)  H_psi(j)=0
 if (pbcx) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanx*k_vecs_r(ik,1))
 if (.not. pbcy  .and. deltany.ne.0)  H_psi(j)=0
 if (pbcy) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltany*k_vecs_r(ik,2))
 if (.not. pbcz .and. deltanz.ne.0)  H_psi(j)=0
 if (pbcz) 	  H_psi(j)=H_psi(j)       *exp(2*pi*(0,1)*deltanz*k_vecs_r(ik,3))

enddo !icf2
enddo	!is2
enddo	!o2

 endif !dist

enddo	!lz
enddo	!ly
enddo	!lx


endif !lreadham



!!!!!!!!!!!!!!!!
!print *, ix,iy,iz, icf, is, H_psi(3*n_nn+2)

!index_psi(:,6)=index_psi(:,1) + index_psi(:,2)*nx + index_psi(:,3)*nx*ny					!site
!index_psi(:,7)=index_psi(:,1) + index_psi(:,2)*nx + index_psi(:,3)*nx*ny + index_psi(:,4)*nx*ny*nz + index_psi(:,5)*2*nx*ny*nz + index_psi(:,8)*4*nx*ny*nz	!site+spin+c/f+o
do i=1,n_max_link
!write (*, "(20i4, i10)"), i,ix,iy,iz,icf,is,o1,index_psi(i,1),index_psi(i,2),index_psi(i,3),index_psi(i,4),index_psi(i,5),index_psi(i,8)
index_psi(i,6)=index_psi(i,1) + index_psi(i,2)*nx + index_psi(i,3)*nx*ny					!site
index_psi(i,7)=index_xyzcso(index_psi(i,1),index_psi(i,2),index_psi(i,3),index_psi(i,4),index_psi(i,5),index_psi(i,8),nx, ny, nz)
enddo


end subroutine ham_psi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_wfc(e_r_tot,w_r_tot)
real (kind=idp),intent(in)		:: e_r_tot(0:total_size-1,0:nk_r-1)
complex(KIND=idpc),intent(in)	 	:: w_r_tot(0:total_size-1,0:total_size-1,0:nk_r-1)
real(KIND=idpc)			:: wx(0:nx-1,3),wy(0:ny-1,3),wz(0:nz-1,3)
integer					::indp1,indp2,indp3,indp4
 character(len=80):: f_e, f_w1, f_w2
 
f_e="(a3,2i6,e15.6)"
f_w1="(3i5,9e13.5)"
f_w2="(i5,9e13.5)"

OPEN(unit=15,file=trim(label)//'/wfc_r',status='unknown')
OPEN(unit=17,file=trim(label)//'/wfcavgx',status='unknown')
OPEN(unit=18,file=trim(label)//'/wfcavgy',status='unknown')
OPEN(unit=19,file=trim(label)//'/wfcavgz',status='unknown')



do ik=0, nk_r-1
do i=0, total_size-1			!eigenvalue
!   write (16, "(i6,2f20.10)") i,e(i)
 if(e_r_tot(i,ik)>e_thr1 .and. e_r_tot(i,ik)<e_thr2) then
  write (15, f_e) "#",ik,i,e_r_tot(i,ik)
 
      do ix=0,nx-1	!x
       do iy=0,ny-1	!y
	 do iz=0,nz-1	!z
    
     indp1=ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz	!c up
     indp2=ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz	!c down
     indp3=ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz	!f up
     indp4=ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz	!f down
		     
     write (15, f_w1) ix,iy,iz,abs(w_r_tot(indp1,i,ik))**2+abs(w_r_tot(indp2,i,ik))**2+abs(w_r_tot(indp3,i,ik))**2+abs(w_r_tot(indp4,i,ik))**2,&
     				 real(w_r_tot(indp1,i,ik)), aimag(w_r_tot(indp1,i,ik)),real(w_r_tot(indp2,i,ik)), aimag(w_r_tot(indp2,i,ik)),&
     				 real(w_r_tot(indp3,i,ik)), aimag(w_r_tot(indp3,i,ik)),real(w_r_tot(indp4,i,ik)), aimag(w_r_tot(indp4,i,ik))
     					

     enddo
    enddo
   enddo
 
  write (15,"(a1)") ""
  write (15,"(a1)") ""

 endif
enddo
enddo

do i=0, total_size-1			!eigenvalue

 
enddo


!!!!!!!!!Average over planes

do ik=0, nk_r-1
do i=0, total_size-1			!eigenvalue

 if(e_r_tot(i,ik)>e_thr1 .and. e_r_tot(i,ik)<e_thr2) then
 
  write (17, f_e) "#",ik, i,e_r_tot(i,ik)
  write (18, f_e) "#",ik, i,e_r_tot(i,ik)
  write (19, f_e) "#",ik, i,e_r_tot(i,ik)
  wx=0
  wy=0
  wz=0

!computes averages


      do ix=0,nx-1	!x
       do iy=0,ny-1	!y
	 do iz=0,nz-1	!z
    
     indp1=ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz	!c up
     indp2=ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz	!c down
     indp3=ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz	!f up
     indp4=ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz	!f down
		     
    			
    					
!     wx(ix,0)=wx(ix,0)+abs(w_r_tot(indp1,i,ik))**2
     wx(ix,1)=wx(ix,1)+abs(w_r_tot(indp2,i,ik))**2
     wx(ix,2)=wx(ix,2)+abs(w_r_tot(indp3,i,ik))**2
     wx(ix,3)=wx(ix,3)+abs(w_r_tot(indp4,i,ik))**2
 !    wy(iy,0)=wy(iy,0)+abs(w_r_tot(indp1,i,ik))**2
     wy(iy,1)=wy(iy,1)+abs(w_r_tot(indp2,i,ik))**2
     wy(iy,2)=wy(iy,2)+abs(w_r_tot(indp3,i,ik))**2
     wy(iy,3)=wy(iy,3)+abs(w_r_tot(indp4,i,ik))**2
  !   wz(iz,0)=wz(iz,0)+abs(w_r_tot(indp1,i,ik))**2
     wz(iz,1)=wz(iz,1)+abs(w_r_tot(indp2,i,ik))**2
     wz(iz,2)=wz(iz,2)+abs(w_r_tot(indp3,i,ik))**2
     wz(iz,3)=wz(iz,3)+abs(w_r_tot(indp4,i,ik))**2
     

     enddo
    enddo
   enddo

!!!writes averages

  do ix=0,nx-1	!x
 !  write (17, f_w2) ix,wx(ix,0)+ wx(ix,1)+ wx(ix,2)+wx(ix,3),wx(ix,0), wx(ix,1), wx(ix,2),wx(ix,3)
  enddo
  do iy=0,ny-1	!x
 !  write (18, f_w2) iy,wy(iy,0)+ wy(iy,1)+ wy(iy,2)+wy(iy,3),wy(iy,0), wy(iy,1), wy(iy,2),wy(iy,3)
  enddo
  do iz=0,nz-1	!x
 !  write (19, f_w2) iz,wz(iz,0)+ wz(iz,1)+ wz(iz,2)+wz(iz,3),wz(iz,0), wz(iz,1), wz(iz,2),wz(iz,3)
  enddo

  write (17,"(a1)") ""
  write (17,"(a1)") ""
  write (18,"(a1)") ""
  write (18,"(a1)") ""
  write (19,"(a1)") ""
  write (19,"(a1)") ""

endif     
enddo
enddo

 close(unit=15)
 close(unit=17)
 close(unit=18)
 close(unit=19)


end subroutine write_wfc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine convert_index(ind, ix,iy,iz,icf,is)
integer	::ind, ix,iy,iz,icf,is

ix=mod(ind,nx)
iy=mod((ind-ix)/nx,ny)
iz=mod((ind-ix-iy*nx)/nx/ny,nz)
icf=mod((ind-ix-iy*nx-iz*nx*ny)/nx/ny/nz,2)
is=mod((ind-ix-iy*nx-iz*nx*ny-icf*nx*ny*nz)/nx/ny/nz/2,2)

!indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
end subroutine convert_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine do_kr_loop
real(kind=idp)	:: time_0,time_1

size_h_kr=2*nz*(l_orbital+l_orbital_f)*nx_kr*ny_kr
!size_h_kr_red=2*nz_scf*(l_orbital+l_orbital_f)*nx_kr*ny_kr
size_h_kr_red=2*nz_g0*(l_orbital+l_orbital_f)*nx_kr*ny_kr
n_k_kr=nkx_kr*nky_kr
n_sites_kr=nz*nx_kr*ny_kr
T=T_min

write(*,"(a25)") ""
write(*,"(a25)") "Kr loop"
write(*,"(a25)") ""

time_0=secnds(0.0)

if (lkred) then
call set_kvecs_kr_pw (k_vecs2d,nk2d,nkx_kr,nky_kr,0,lkshifted)	!reduced set
else
call set_kvecs2d (k_vecs2d,nk2d, nkx_kr, nky_kr,.not. lkshifted)				!whole set
endif

!k_vecs2d(:,1)=k_vecs2d(:,1)/nx_kr
!k_vecs2d(:,2)=k_vecs2d(:,2)/ny_kr

 call allocate_var_kr

if (.not. lreadham) then

! call set_phi_k(phi_kr_sigma_alpha_tot,phi_kr_alpha_sigma_tot, kin_energy_kr_tot, delta_kr_tot,nk2d,k_vecs2d,2)
! call set_kin_energy2_kr(kin_energy2_kr_tot,phi_kr_sigma_alpha2_tot,phi_kr_alpha_sigma2_tot,nk2d,k_vecs2d)
! call set_kin_energy3_kr(kin_energy3_kr_tot,phi_kr_sigma_alpha3_tot,phi_kr_alpha_sigma3_tot,nk2d,k_vecs2d)

endif

!call factor_k(factor_kr_tot,nk2d,k_vecs2d,2)

zero_kr=0.0d0
one_kr=1.0d0

if (lscfkr) call build_v_kr(v_kr_tot)
if (lscfkr) call build_V_kr2(h_kr_psi_tot, index_kr_psi_tot)


if (lscfkr) call write_k_points_kr(k_vecs2d)	

 call set_initialvalues_kr

if (lscfkr) then

print *, "SCF kr loop"

time_1=secnds(0.0)
time_kr_setvar=time_1-time_0


!scf-loop
 diff=1
 i_iter=0
 ok=.false.

 do while (diff>threshold_kr .and. lscfloop_kr)
  
  i_iter=i_iter+1
  print *, i_iter


  e_kr_tot=0		!vector of all eigenvalues
  w_kr_tot=0
 
! k-vecs loop
  do ik=0,nk2d-1
   time_0=secnds(0.0)
 
 !  if (i_iter==1)  write(*,"(i5,a15,2(a10,f10.5))") ik, "k point",  "Delta_k", delta_kr_tot(ik)
  
!if (.not. lreadham)   call build_ham_kr(ik,Ham_kr,b_kr,lambda_kr, mu_kr,.true.,nk2d,phi_kr_sigma_alpha_tot,phi_kr_alpha_sigma_tot, kin_energy_kr_tot, &
!   			kin_energy2_kr_tot, kin_energy3_kr_tot,phi_kr_sigma_alpha2_tot, phi_kr_alpha_sigma2_tot,phi_kr_sigma_alpha3_tot, phi_kr_alpha_sigma3_tot, k_vecs2d)
!if (lreadham)         call build_ham_kr_read(ik,Ham_kr,b_kr,lambda_kr, mu_kr,.true.,nk2d,k_vecs2d)
		      call build_ham_kr2(ik,Ham_kr,b_kr, lambda_kr,mu_kr,k_vecs2d,nk2d, .true.)



!   if (n_sites_kr< 10 .and. i_iter==1)     call write_Ham_kr(ham_kr)
 
 !   call build_ham_kr(ik,V_kr,1.d0,0.d0, 0.d0)
 !  do i=0,3
 ! V_k(i,i)=0
  !enddo 


!   if (n_k< 5 .and. i_iter==1)     call write_Ham_k(ham_k)

  time_1=secnds(0.0)
  time_kr_buildh=time_kr_buildh+time_1-time_0

  time_0=secnds(0.0)

   call diagonalize_ham_kr(Ham_kr,e_kr,w_kr)


!  if (n_k< 5 .and. i_iter==1)    call write_energies_k (e_k,w_k)

   e_kr_tot(:,ik)=e_kr(:)
   w_kr_tot(:,:,ik)= w_kr(:,:)

  time_1=secnds(0.0)
  time_kr_diag=time_kr_diag+time_1-time_0

  enddo		!k points

!!!!!!!!!!!!Compute chemical potential mu

  time_0=secnds(0.0)

  if (lfixmu_kr) then
   mu_new_kr=mu_kr
  else  
   call compute_mu_kr(e_kr_tot)
!   call compute_mu_kr_1(e_kr_tot)
  endif
  
  time_1=secnds(0.0)
  time_kr_mu=time_kr_mu+time_1-time_0

!!!!!!!!Compute occupations, new b/lambda, and total free energy

  time_0=secnds(0.0)

!only 1x1 cell, faster
if (nx_kr*ny_kr==1)    call compute_nf_kr (nf_kr,nc_kr,b_new_kr,e_kr_tot,w_kr_tot,ef_kr, lambda_new_kr, nfc_kr,diff,nff_kr,0.0d0)		!ottengo b_new old method; ok
!larger cell, slower for 1x1 cell
if (nx_kr*ny_kr .ne. 1)&
 call compute_nf_r  (nf_kr,nc_kr,nfc_kr,nff_kr,b_new_kr,lambda_new_kr,b_kr,lambda_kr,0.0d0,e_kr_tot,w_kr_tot,diff,k_vecs2d,nk2d,n_k_kr, size_h_kr,n_sites_kr,nx_kr,ny_kr,nz,h_kr_psi_tot,index_kr_psi_tot,lvac_kr, lambda_par_kr)


!writes to file, not to be read, but possibly substituted to b_kr2 to restart a SCF calculation
  OPEN(unit=131,file=trim(label)//'/b_kr3',status='unknown')
 write(131,"(3i4)") nx_kr,ny_kr,nz

 do ix=0,nx_kr-1	!x
   do iy=0,ny_kr-1	!y
    do iz=0,nz-1	!z
 ind=index_xyz(ix,iy,iz,nx_kr,ny_kr)
 write(131,"(3i4,100f15.10)") ix,iy,iz, b_kr(ind), lambda_kr(ind)
   enddo
  enddo
 enddo
 write(131,"(f20.10)")  mu_kr
 close (unit=31)



!!!!!!!!!!!Check selfconsistency

!diff=0

!do i=0, nz-1
!  diff=diff+(1-b_new_kr(i)**2-nf_kr(i))**2/nz
!enddo

!if (.not. lfixmu) diff=diff+((Nfctot-N_el_kr)/n_k_kr)**2

!diff=sqrt(diff)


!!!!!!!!!!!!!Set new b and lambda

  if (diff<threshold_kr) then
  ok=.true.
  else
  b_kr=alpha_kr*b_new_kr+(1-alpha_kr)*b_kr
  lambda_kr=alpha_kr*lambda_new_kr+(1-alpha_kr)*lambda_kr
  mu_kr=alpha_kr*mu_new_kr+(1-alpha_kr)*mu_kr
  endif

  time_1=secnds(0.0)
  time_kr_nf=time_kr_nf+time_1-time_0

 enddo	!SCF loop



! nkr_plot=(nkx_plot+1)*(nky_plot+1)

endif	!lscfkr



 time_0=secnds(0.0)

 if (ok) then
ef_kr=0

 write(*,"(a20)")  " "
write (22, "(a5,21(a15))"), "#i", "b","bnfc" ,"nc","nf", "ef",  "lambda",  "nfc","nff","nfk","nfck"
do i =0, nz-1
write(22,"(i5,100f15.10)")  i, b_kr(i), (nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(2*lambda_kr(i)-2*nfk_kr(i)),nc_kr(i),nf_kr(i), ef_kr(i), lambda_kr(i),  nfc_kr(i), nff_kr(i),nfk_kr(i),nfck_kr(i)!, E_tot_kr
enddo
 write(*,"(a20)")  " "
write(22,"(a10,i10)")  "#n_iter=", i_iter

 if (lscfkr) then
  OPEN(unit=31,file=trim(label)//'/b_kr2',status='unknown')
 write(31,"(3i4)") nx_kr,ny_kr,nz

 do ix=0,nx_kr-1	!x
   do iy=0,ny_kr-1	!y
    do iz=0,nz-1	!z
 ind=index_xyz(ix,iy,iz,nx_kr,ny_kr)
 write(31,"(3i4,100f15.10)") ix,iy,iz, b_kr(ind), lambda_kr(ind)
   enddo
  enddo
 enddo
 write(31,"(f20.10)")  mu_kr
 close (unit=31)
 endif


!print *, nff_kr




 if (loopk) write(*,"(a20,6f14.8,i5)") "K loop: ", T,b_new_k, lambda_new_k, mu_new_k, E_tot_k


if (.not. lscfkr) mu_new_kr=mu_kr
e_kr_tot=e_kr_tot-mu_new_kr

if (lscfkr) call write_bands_kr(e_kr_tot)

!if (lwfc) call write_wfc_kr(e_kr_tot,w_kr_tot)

!if (lspectrfunct_kr) call spectralfunction_kr(e_kr_tot,w_kr_tot)

 endif

 time_0=secnds(0.0)

 if (lnscfk) call set_kr_vecs_nscf(kr_vecs_nscf, nkx_plot, nky_plot)
 if (lnscfk) call do_kr_nscf(b_kr,lambda_kr, mu_kr)

 time_1=secnds(0.0)
 time_nscfkr=time_1-time_0


time_1=secnds(0.0)
time_kr_write=time_1-time_0



!call deallocate_var_kr

write(*,"(a25)") ""
write(*,"(a25)") "end kr loop"
write(*,"(a25)") ""



end subroutine do_kr_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_var_kr

!allocate(phi_k_sigma_alpha_tot(0:n_k-1,0:1,0:1))
!allocate(phi_k_alpha_sigma_tot(0:n_k-1,0:1,0:1))
if (.not. lreadham) allocate(phi_kr_sigma_alpha_tot(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1))
if (.not. lreadham) allocate(phi_kr_alpha_sigma_tot(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk2d-1))
if (.not. lreadham) allocate(phi_kr_sigma_alpha2_tot(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1,n_nn2))
if (.not. lreadham) allocate(phi_kr_alpha_sigma2_tot(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk2d-1,n_nn2))
if (.not. lreadham) allocate(phi_kr_sigma_alpha3_tot(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1,n_nn3))
if (.not. lreadham) allocate(phi_kr_alpha_sigma3_tot(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk2d-1,n_nn3))
if (.not. lreadham) allocate(kin_energy_kr_tot(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1,0:1))
if (.not. lreadham) allocate(kin_energy2_kr_tot(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1,n_nn2,0:1))
if (.not. lreadham) allocate(kin_energy3_kr_tot(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk2d-1,n_nn3,0:1))
allocate(ham_kr(0:size_h_kr-1,0:size_h_kr-1))
allocate(e_kr(0:size_h_kr-1))
allocate(e_kr_tot(0:size_h_kr-1,0:nk2d-1))
allocate(w_kr(0:size_h_kr-1,0:size_h_kr-1))
allocate(w_kr_rec(0:size_h_kr-1,0:size_h_kr-1))
allocate(w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1))
allocate(V_kr(0:size_h_kr-1,0:size_h_kr-1))
allocate(v_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1))
allocate(k_vecs_kr(0:n_k_kr-1,1:4))
allocate(ef_kr(0:n_sites_kr-1))
allocate(b_kr(0:n_sites_kr-1))
allocate(one_kr(0:n_sites_kr-1))
allocate(zero_kr(0:n_sites_kr-1))
allocate(b_new_kr(0:n_sites_kr-1))
allocate(lambda_kr(0:n_sites_kr-1))
allocate(lambda_new_kr(0:n_sites_kr-1))
allocate(nf_kr(0:n_sites_kr-1))	
allocate(nfk_kr(0:n_sites_kr-1))	
allocate(nc_kr(0:n_sites_kr-1))
allocate(nfc_kr(0:n_sites_kr-1))
allocate(nfck_kr(0:n_sites_kr-1))
allocate(nff_kr(0:n_sites_kr-1))
allocate(delta_kr_tot(0:nk2d-1))
!allocate(factor_kr_tot(0:nk2d-1))
allocate(lvac_kr(0:n_sites_kr-1))
allocate(ef_site_kr(0:n_sites_kr-1,0:l_orbital_f-1))
allocate(ec_site_kr(0:n_sites_kr-1,0:max(l_orbital,1)-1))
allocate(index_kr_psi_tot(0:size_h_kr-1,n_max_link,3))


!allocate(index_psi(n_max_link,7))
!allocate(h_psi(n_max_link))
!allocate(index_psi_tot(0:total_size-1,n_max_link,3))
!allocate(h_psi_tot(0:total_size-1,n_max_link))

end subroutine allocate_var_kr

subroutine deallocate_var_kr

if (.not. lreadham) deallocate(phi_kr_sigma_alpha_tot)
if (.not. lreadham) deallocate(phi_kr_alpha_sigma_tot)
if (.not. lreadham) deallocate(phi_kr_sigma_alpha2_tot)
if (.not. lreadham) deallocate(phi_kr_alpha_sigma2_tot)
if (.not. lreadham) deallocate(phi_kr_sigma_alpha3_tot)
if (.not. lreadham) deallocate(phi_kr_alpha_sigma3_tot)
if (.not. lreadham) deallocate(kin_energy_kr_tot)
if (.not. lreadham) deallocate(kin_energy2_kr_tot)
if (.not. lreadham) deallocate(kin_energy3_kr_tot)
deallocate(ham_kr)
deallocate(e_kr)
deallocate(e_kr_tot)
deallocate(w_kr)
deallocate(w_kr_tot)
deallocate(V_kr)
deallocate(k_vecs_kr)
deallocate(ef_kr)
deallocate(b_kr)
deallocate(one_kr)
deallocate(zero_kr)
!deallocate(b_new_kr)
deallocate(lambda_kr)
!deallocate(lambda_new_kr)
deallocate(nf_kr)	
deallocate(nfk_kr)	
deallocate(nc_kr)
deallocate(nfc_kr)
deallocate(nfck_kr)
deallocate(nff_kr)
deallocate(v_kr_tot)
deallocate(lvac_kr)
deallocate(ef_site_kr)
deallocate(ec_site_kr)
deallocate(H_kr_psi_tot)
deallocate(index_kr_psi_tot)

end subroutine deallocate_var_kr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine convert_index_k_kr(ind, ikx,iky)
integer,intent(out)	::ikx,iky
integer,intent(in)	::ind


ikx=mod(ind,nkx)
iky=mod((ind-ikx)/nkx,nky)
!ikz=mod((ind-ikx-iky*nkx)/nkx/nky,nkz)

!ind=ix + iy*nx + iz*nx*ny



end subroutine convert_index_k_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_k_points_kr(k_vecs2d)
real(kind=idp), intent(in)		::k_vecs2d(0:nk2d-1,1:4)
integer	::i

do i=6,24,18

 write(i,"(a15)") ""

 write(i,"(a20)")  "List of K_(x,y)-points:"
 write(i,"(5a10)")  "kx", "ky",  "F_k", "Delta_k^2"
 write(i,"(a35)") "Phi_sigma_alpha, Phi_alpha_sigma"
 write(i,"(a15)") ""

 do ik=0, nk2d-1 

  write(i,"(i8,6f10.5)")  ik, 	 k_vecs2d(ik,1), k_vecs2d(ik,2), k_vecs2d(ik,3),k_vecs2d(ik,4) !, delta_kr_tot(ik)

!if (.not. lreadham) then
!     do is=0,1
!     write(i,"(4 ('(',f6.2,',',f6.2,')','  ') )")  real(phi_kr_sigma_alpha_tot(0,0,is,0,ik)),aimag(phi_kr_sigma_alpha_tot(0,0,is,0,ik)), &
!     						   real(phi_kr_sigma_alpha_tot(0,0,is,1,ik)),aimag(phi_kr_sigma_alpha_tot(0,0,is,1,ik)),&
!   						   real(phi_kr_alpha_sigma_tot(0,0,is,0,ik)),aimag(phi_kr_alpha_sigma_tot(0,0,is,0,ik)), &
!   						   real(phi_kr_alpha_sigma_tot(0,0,is,1,ik)),aimag(phi_kr_alpha_sigma_tot(0,0,is,1,ik))
!    enddo 
!endif
!  write(i, " ")
 enddo


 write(i,"(a10)") "**********"

enddo

end subroutine write_k_points_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sum_i(nn) exp(i*k*r_i)
subroutine factor_Kr(ik,f_kr)
integer, intent(in)	::ik
real (kind=idp),intent(out)	::f_kr
integer			::i, n_nn_xy
complex(kind=idpc)	::f_k_c

f_k_c=0
if(lkz) n_nn_xy=n_nn-2
if(.not. lkz) n_nn_xy=n_nn

do i=1, n_nn_xy
f_k_c=f_k_c+exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs_kr(ik,1)+nn_ijk(i,2)*k_vecs_kr(ik,2)))
!print *, nn_ijk(i,1), k_vecs(ik,1), nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3),&
!(2*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))!,exp(2*(0,1)*pi*(nn_ijk(i,1)*k_vecs(ik,1)+nn_ijk(i,2)*k_vecs(ik,2)+nn_ijk(i,3)*k_vecs(ik,3)))
enddo


if (aimag(f_k_c)>1e-9) print *, "warning f_c_k", f_k_c

f_kr=real(f_k_c)

end subroutine factor_Kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kr(ik,Ham_kr,b_kr,lambda_kr, mu_kr, ldiag,nk,phi_sa, phi_as, ke, ke2, ke3, phi_sa2, phi_as2, phi_sa3, phi_as3, k_vecs)
complex(KIND=idpc),intent(out) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1)
complex(KIND=idpc),intent(in) 	:: phi_sa(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1),phi_as(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1)
complex(KIND=idpc),intent(in) 	:: phi_sa2(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn2),phi_as2(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1,n_nn2)
complex(KIND=idpc),intent(in) 	:: phi_sa3(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn3),phi_as3(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk-1,n_nn3)
complex(KIND=idpc),intent(in) 	:: ke(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,0:1),ke2(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn2,0:1),&
				ke3(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk-1,n_nn3,0:1)
real(KIND=idp),intent(in) 	:: b_kr(0:n_sites_kr-1), lambda_kr(0:n_sites_kr-1), mu_kr,k_vecs(0:nk-1,6)
integer, intent(in)		::ik,nk
integer				:: i
logical, intent(in)		::ldiag


Ham_kr=0

ix=0
iy=0
do iz=0, nz-1
! do ix=0, nx_kr-1
! do iy=0, ny_kr-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
    do icf2=0,1
     do o2=0,l_orb(icf2)-1
      do is2=0,1
       do iz2=0, nz-1
! 	do ix2=0, nx_kr-1
!	 do iy2=0, ny_kr-1
 
 indp=index_xyzcso (ix, iy, iz, icf, is, o1, nx_kr, ny_kr, nz)
 ind =index_xyz(ix,iy,iz,nx_kr, ny_kr)

 ind2p=index_xyzcso (ix2, iy2, iz2, icf2, is2, o2, nx_kr, ny_kr, nz)
 ind2 =index_xyz(ix2,iy2,iz2,nx_kr, ny_kr)


 
!onsite energy
if (indp==ind2p .and. icf==0 .and. ldiag) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+ec_site_kr(ind,o1)-mu_kr					!c
if (indp==ind2p .and. icf==1 .and. ldiag) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+ef_site_kr(ind,o1)-lambda_kr(ind)-mu_kr		!f

if (iz==iz2) then
!ke
if (icf==0 .and. icf2==0) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+tc_*ke(o1,o2,is,is2,ik,0)		!contains 2nd nn in the plane and 3rd
if (icf==1 .and. icf2==1) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+b_kr(ind)**2*tf_*ke(o1,o2,is,is2,ik,1)  !contains 2nd nn in the plane and 3rd

!sigma_alpha
if (icf==0 .and. icf2==1) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+phi_sa(o1,o2,is,is2,ik)*Vcf_*b_kr(ind)	!contains 2nd nn in the plane and 3rd
!alpha_sigma
if (icf==1 .and. icf2==0) Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+phi_as(o1,o2,is,is2,ik)*Vcf_*b_kr(ind)	!contains 2nd nn in the plane and 3rd
endif
!onsite hybr
if (icf==0 .and. icf2==1 .and. phi_type=="ons") Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+Vcf_*b_kr(ind)					!c
if (icf==1 .and. icf2==0 .and. phi_type=="ons") Ham_kr(indp,ind2p)=Ham_kr(indp,ind2p)+Vcf_*b_kr(ind)		!f



!if (icf==0 .and. icf2==1 .and. iz2==iz) Ham_kr(indp,ind2p)=Vcf_*phi_sa(is,is2,ik)*b_kr(ind)
!if (icf==1 .and. icf2==0 .and. iz2==iz) Ham_kr(indp,ind2p)=Vcf_*phi_as(is,is2,ik)*b_kr(ind)

!between planes
!if there are at least 3 sites along z: check for 2 sites!
!1st nn
      do i=n_nn-1,n_nn			!only along z
       if ((iz==iz2+nn_ijk(i,3).and. lkz) .or. &
       ((pbcz .and. nz>2 .and. iz==iz2+nn_ijk(i,3)-nz .and. lkz)) .or. &
       (pbcz .and. nz>2 .and. iz==iz2+nn_ijk(i,3)+nz .and. lkz)) then	
        if (icf==0 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tc_*kin_energy(o1,o2,is,is2,i,0)
        if (icf==1 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tf_*kin_energy(o1,o2,is,is2,i,1)*b_kr(ind)*b_kr(ind2)
        if (icf==0 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_sigma_alpha(o1,o2,is,is2,i)*b_kr(ind2)
        if (icf==1 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_alpha_sigma(o1,o2,is,is2,i)*b_kr(ind)
!        if (icf==0 .and. icf2==0 .and. is==is2) Ham_kr(indp, ind2p)=tc_ 
!        if (icf==1 .and. icf2==1 .and. is==is2) Ham_kr(indp, ind2p)=tf_*b_kr(ind)*b_kr(ind2)
!        if (icf==0 .and. icf2==1 .and. Vtki)    Ham_kr(indp, ind2p)=Vcf_*phi_sigma_alpha(0,0,is,is2,i)*b_kr(ind2)
!        if (icf==1 .and. icf2==0 .and. Vtki)    Ham_kr(indp, ind2p)=Vcf_*phi_alpha_sigma(0,0,is,is2,i)*b_kr(ind)
       endif
      
!       if (pbcz .and. nz>2 .and. iz2==iz+nn_ijk(i,3)-nz .and. lkz)  then
!        if (icf==0 .and. icf2==0 .and. is==is2) Ham_kr(indp, ind2p)=tc_ 
!        if (icf==1 .and. icf2==1 .and. is==is2) Ham_kr(indp, ind2p)=tf_*b_kr(ind)*b_kr(ind2)
!        if (icf==0 .and. icf2==1 .and. Vtki) Ham_kr(indp, ind2p)=Vcf_*phi_sigma_alpha(0,0,is,is2,i)*b_kr(ind2)
!        if (icf==1 .and. icf2==0 .and. Vtki) Ham_kr(indp, ind2p)=Vcf_*phi_alpha_sigma(0,0,is,is2,i)*b_kr(ind)
!       endif

!       if (pbcz .and. nz>2 .and. iz2==iz+nn_ijk(i,3)+nz .and. lkz)  then
!        if (icf==0 .and. icf2==0 .and. is==is2) Ham_kr(indp, ind2p)=tc_ 
!        if (icf==1 .and. icf2==1 .and. is==is2) Ham_kr(indp, ind2p)=tf_*b_kr(ind)*b_kr(ind2)
!        if (icf==0 .and. icf2==1 .and. Vtki) Ham_kr(indp, ind2p)=Vcf_*phi_sigma_alpha(0,0,is,is2,i)*b_kr(ind2)
!        if (icf==1 .and. icf2==0 .and. Vtki) Ham_kr(indp, ind2p)=Vcf_*phi_alpha_sigma(0,0,is,is2,i)*b_kr(ind)
!       endif
      enddo
     
!2nd nn
      do i=1,n_nn2			
      if(nn2_ijk(i,3)==1 .or. nn2_ijk(i,3)==-1) then !only along z 
       if ((iz==iz2+nn2_ijk(i,3).and. lkz2) .or. &
       ((pbcz .and. nz>2 .and. iz==iz2+nn2_ijk(i,3)-nz .and. lkz2)) .or. &
       (pbcz .and. nz>2 .and. iz==iz2+nn2_ijk(i,3)+nz .and. lkz2)) then	
        if (icf==0 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tc_*ke2(o1,o2,is,is2,ik,i,0)
        if (icf==1 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tf_*ke2(o1,o2,is,is2,ik,i,1)*b_kr(ind)*b_kr(ind2)
        if (icf==0 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_sa2(o1,o2,is,is2,ik,i)*b_kr(ind2)
        if (icf==1 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_as2(o1,o2,is,is2,ik,i)*b_kr(ind)
       endif
       endif
      enddo
!3rd nn
      do i=1,n_nn3			
      if(nn3_ijk(i,3)==1 .or. nn3_ijk(i,3)==-1) then !only along z 
       if ((iz==iz2+nn3_ijk(i,3).and. lkz2) .or. &
       ((pbcz .and. nz>2 .and. iz==iz2+nn3_ijk(i,3)-nz .and. lkz2)) .or. &
       (pbcz .and. nz>2 .and. iz==iz2+nn3_ijk(i,3)+nz .and. lkz2)) then	
        if (icf==0 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tc_*ke3(o1,o2,is,is2,ik,i,0)
        if (icf==1 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+tf_*ke3(o1,o2,is,is2,ik,i,1)*b_kr(ind)*b_kr(ind2)
        if (icf==0 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_sa3(o1,o2,is,is2,ik,i)*b_kr(ind2)
        if (icf==1 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Vcf_*phi_as3(o1,o2,is,is2,ik,i)*b_kr(ind)
       endif
       endif
      enddo

     
     
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
enddo
enddo

!print *, ham_kr

!sets diagonal to zero
!if (.not. ldiag) then
!do i=0,size_h_kr-1
!Ham_kr(i,i)=0
!enddo
!endif

!Ham(ri,r'j)=<ri|H|r'j>=<r-r'i|H|j>==<l'i|H|j>--> r=r'+l



!if (n_sites <5) then
do i=0, size_h_kr-1
do j=i, size_h_kr-1
if (abs(Ham_kr(j,i)-conjg(Ham_kr(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham_kr(i,j), Ham_kr(j,i)
stop
endif
enddo
enddo
!endif

end subroutine build_ham_kr



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kr2(ik,Ham,b, lambda,mu,k_vecs,nk,ldiag)
complex(KIND=idpc),intent(out) 	:: Ham(0:size_h_kr-1,0:size_h_kr-1)
real(KIND=idp),intent(in) 	:: b(0:n_sites_kr-1), lambda(0:n_sites_kr-1), mu,k_vecs(0:nk-1,4)
integer, intent(in)::ik,nk
logical, intent(in):: ldiag

ham=0

do ix=0,nx_kr-1	!x
 do iy=0,ny_kr-1	!y
  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1
     
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx_kr, ny_kr, nz)
	index_psi=0
	h_psi=0
      call ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,b,lambda,mu,ldiag,nx_kr,ny_kr,nz,ec_site_kr,ef_site_kr, k_vecs,nk)
!	   ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,b,lambda,mu,ldiag,nx,   ny,   nz,ec_site   ,ef_site,    k_vecs_r, n_k)

     do i=1,n_max_link
      ind2p=index_psi(i,7)
      Ham(indp,ind2p)=Ham(indp,ind2p)+h_psi(i)
!write (*, "(2i4,2f10.5)") i,ind2p,real(h_psi(i)), aimag(h_psi(i))
     enddo

     enddo
    enddo
   enddo
  enddo
 enddo
enddo




do i=0, size_h_kr-1
do j=i, size_h_kr-1
if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham(i,j), Ham(j,i)
stop
endif
enddo
enddo

end subroutine build_ham_kr2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kr_kdotp(ik,Ham,b, lambda,mu,k_vecs,nk,ldiag)
complex(KIND=idpc),intent(out) 	:: Ham(0:size_h_kr-1,0:size_h_kr-1)
real(KIND=idp),intent(in) 	:: b(0:n_sites_kr-1), lambda(0:n_sites_kr-1), mu,k_vecs(0:nk-1,6)
integer, intent(in)::ik,nk
logical, intent(in):: ldiag
real(KIND=idp)::sincos,sincos2, disp_n
complex(KIND=idpc)::kp, km


!k_z -> -i d/dz

!print *, "build ham kr kdp"
ham=0

kp=(k_vecs(ik,1)+(0,1)*k_vecs(ik,2))*Vcf_*b_k*2*pi	!kx+iky
km=(k_vecs(ik,1)-(0,1)*k_vecs(ik,2))*Vcf_*b_k*2*pi	!kx-iky
kpar2=(k_vecs(ik,1)**2+k_vecs(ik,2)**2)*4*pi**2		!kx^2+ky^2
k2xy=(-k_vecs(ik,1)**2+k_vecs(ik,2)**2)*4*pi**2		!ky^2-kx^2


!ham_k(0,0)=eckp1-mu_k+tc_*kz2*g1d+tc_*kpar2*l1d	!e_d1
!ham_k(1,1)=ham_k(0,0)
!ham_k(2,2)=eckp2-mu_k+tc_*kz2*g2d+tc_*kpar2*l2d	!e_d2
!ham_k(3,3)=ham_k(2,2)
!ham_k(4,4)=efkp1-mu_k-lambda_k+tf_*b_k**2*kz2*g1f+tf_*b_k**2*kpar2*l1f	!ef1
!ham_k(5,5)=ham_k(4,4)
!ham_k(6,6)=efkp2-mu_k-lambda_k+tf_*b_k**2*kz2*g2f+tf_*b_k**2*kpar2*l2f	!e_f2
!ham_k(7,7)=ham_k(6,6)
!if (l_orbital_f==3) ham_k(8,8)=efkp7-mu_k-lambda_k+tf_*b_k**2*kz2*g7+tf_*b_k**2*l7*kpar2	!e_f7
!if (l_orbital_f==3) ham_k(9,9)=ham_k(8,8)
!non diag kin terms
!ham_k(0,2)=tc_*td12*k2xy
!ham_k(1,3)=ham_k(0,2)
!ham_k(4,6)=tf_*b_k**2*tf12*k2xy
!ham_k(5,7)=ham_k(4,6)


!here iz denotes n, the index of the wfc sin [pi n z/(N+1)], n,z=1,..,N


do ix=0,nx_kr-1	!x
 do iy=0,ny_kr-1	!y
  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1
          
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx_kr, ny_kr, nz)

!continouous model does not work
!disp_n=(iz+1)**2/(nz+1)**2*pi**2	!free dispersion
!discretized model
if (.not. pbcz) disp_n=2-2*cos(pi*(iz+1)/(nz+1))	!discretized dispersion
if (      pbcz) disp_n=2-2*cos(pi*(iz)/(nz))	!discretized dispersion


!onsite enrgy, kx^2+ky^2, dz^2
!iz k index of sin(pi*x*k)/(N+1)
if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*kpar2*l1d+tc_*g1d*disp_n
if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*kpar2*l2d+tc_*g2d*disp_n
if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*kpar2*l1f+tf_*b_k**2*g1f*disp_n
if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*kpar2*l2f+tf_*b_k**2*g2f*disp_n
if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tf_*b_k**2*kpar2*l7 +tf_*b_k**2*g7*disp_n
!iz layer index
!if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*kpar2*l1d+2*tc_*g1d
!if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*kpar2*l2d+2*tc_*g2d
!if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*kpar2*l1f+2*tf_*b_k**2*g1f
!if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*kpar2*l2f+2*tf_*b_k**2*g2f
!if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tc_*b_k**2*kpar2*l7+ 2*tf_*b_k**2*g7
!m78


do ix2=0,nx_kr-1	!x
 do iy2=0,ny_kr-1	!y
  do iz2=0,nz-1	!z
   do icf2=0,1	!c/f
    do is2=0,1	!up/down 
     do o2=0,l_orb(icf2)-1
      ind2p=index_xyzcso (ix2 , iy2 , iz2, icf2, is2, o2, nx_kr, ny_kr, nz)

!non diag ke (perturbation)
if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy	!td12
if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy	!td12
if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78
if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78


!kinetic energy for discretized model --> iz is layer index
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tc_*g1d	!d1
!if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tc_*g2d	!d2
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g1f	!f1
!if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g2f	!f2
!if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==2 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g7	!f7

!hybridization for discretized model --> iz is layer index
!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=+f1*(2*is-1)*Vcf_*b_k/2!*pi		!df1
!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=-f1*(2*is-1)*Vcf_*b_k/2!*pi		!df1
!if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=+f2*(2*is-1)*Vcf_*b_k/2!*pi		!df2
!if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=-f2*(2*is-1)*Vcf_*b_k/2!*pi		!df2
!if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=-f1*(2*is-1)*Vcf_*b_k/2!*pi		!fd1
!if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=+f1*(2*is-1)*Vcf_*b_k/2!*pi		!fd1
!if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=-f2*(2*is-1)*Vcf_*b_k/2!*pi		!fd2
!if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=+f2*(2*is-1)*Vcf_*b_k/2!*pi		!fd2



!if (iz.ne.iz2)sincos= 2*(-1+(-1)**(iz+iz2))*(iz2+1)*(iz+1)/((iz+1)**2-(iz2+1)**2)/(nz+1)
!if (iz.ne.iz2)sincos= 2*(-1+(-1)**(iz+iz2))*sin(pi*(iz2+1)*(iz+1)/(nz+1))/((iz+1)**2-(iz2+1)**2)/pi
!if (iz.ne.iz2)sincos=2*sin(pi*(iz+1)*(iz2+1)/(nz+1))	!discretized hybr (try)

!continouos solution does not work
!if (iz.ne.iz2) sincos= 2*(-1+(-1)**(iz+iz2))*(iz2+1)*(iz+1)/((iz+1)**2-(iz2+1)**2)/(nz+1)
!discretized solution
if (.not. pbcz) sincos=hybr_k(iz2+1,iz+1,nz)*norm_k(iz+1,nz)*norm_k(iz2+1,nz)
if (      pbcz) sincos=hybr_k_pbc(iz2,iz,nz)
!if (      pbcz .and. iz==iz2) sincos=-(0,1)*sin(2*pi*k/nz)  !sincos=hybr_k_pbc(iz2,iz,nz)


!hybr kz-> -i dz (belongs to H_0)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=-f1*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=-f2*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!df2
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=-f7*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!df7
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=f1*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!fd1
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=f2*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!fd2
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2.ne.iz) ham(indp,ind2p)=f7*sincos*(2*is-1)*Vcf_*b_k!*2*pi		!fd7


!if (iz.ne.iz2)sincos2=2*(-1+(-1)**(iz2+iz))*(iz+1)*(iz2+1)/((iz2+1)**2-(iz+1)**2)/(nz+1)
!if (iz.ne.iz2)sincos2=-sincos

!kin energy t12 k^2_y-k^2_x, no dz term
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf2==0 .and. icf==0 .and. o2==0 .and. o1==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf2==0 .and. icf==0 .and. o2==1 .and. o1==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf2==1 .and. icf==1 .and. o2==0 .and. o1==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf2==1 .and. icf==1 .and. o2==1 .and. o1==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!hybr k+, k- (perturbation)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h1*km	!df
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h1*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*kp
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h2*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h2*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h7*km	!df7
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h7*kp
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h72*kp	!df72
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h72*km

if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h1*kp	!fd
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h1*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*km
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h2*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h2*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h7*kp	!fd7
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h7*km
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h72*km	!fd72
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h72*kp
!add gamma_7

     enddo
    enddo
   enddo
  enddo
 enddo
enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo




do i=0, size_h_kr-1
do j=i, size_h_kr-1
!if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9 .and. abs(Ham(j,i)+conjg(Ham(i,j)))>1e-9) then
if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham(i,j), Ham(j,i)
!stop
endif
enddo
enddo

end subroutine build_ham_kr_kdotp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kr_kdotp_x(ik,Ham,b, lambda,mu,k_vecs,nk,ldiag)
complex(KIND=idpc),intent(out) 	:: Ham(0:size_h_kr-1,0:size_h_kr-1)
real(KIND=idp),intent(in) 	:: b(0:n_sites_kr-1), lambda(0:n_sites_kr-1), mu,k_vecs(0:nk-1,6)
integer, intent(in)::ik,nk
logical, intent(in):: ldiag
real(KIND=idp)::sincos,sincos2, disp_n,ky2,kz2,kz
complex(KIND=idpc)::kp, km

!k_x -> -i d/dx

!print *, "build ham kr kdp"
ham=0


!here kx <-->kz
kp=+(0,1)*k_vecs(ik,2)*Vcf_*b_k*2*pi	!+iky
km=-(0,1)*k_vecs(ik,2)*Vcf_*b_k*2*pi	!-iky
kz=(k_vecs(ik,1)-0.5)*Vcf_*b_k*2*pi	!kz
ky2=k_vecs(ik,2)**2*4*pi**2		!ky^2
kz2=(k_vecs(ik,1)-0.5)**2*4*pi**2	!kz^2


!here iz denotes n, the index of the wfc sin [pi n z/(N+1)], n,z=1,..,N


do ix=0,nx_kr-1	!x
 do iy=0,ny_kr-1	!y
  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1
          
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx_kr, ny_kr, nz)

!disp_n=(iz+1)**2/(nz+1)**2*pi**2	!free dispersion
disp_n=2-2*cos(pi*(iz+1)/(nz+1))	!discretized dispersion

!onsite enrgy, kx^2+ky^2, dz^2
!iz k index of sin(pi*x*k)/(N+1)
if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*(ky2*l1d+kz2*g1d)+tc_*l1d*disp_n
if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*(ky2*l2d+kz2*g2d)+tc_*l2d*disp_n
if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*(ky2*l1f+kz2*g1f)+tf_*b_k**2*l1f*disp_n
if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*(ky2*l2f+kz2*g2f)+tf_*b_k**2*l2f*disp_n
if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tf_*b_k**2*(ky2*l7+kz2*g7)  +tf_*b_k**2*l7*disp_n
!if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tf_*b_k**2*kpar2*l7 +tf_*b_k**2*g7*disp_n
!iz layer index
!if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*kpar2*l1d+2*tc_*g1d
!if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*kpar2*l2d+2*tc_*g2d
!if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*kpar2*l1f+2*tf_*b_k**2*g1f
!if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*kpar2*l2f+2*tf_*b_k**2*g2f
!if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tc_*b_k**2*kpar2*l7+ 2*tf_*b_k**2*g7
!m78


do ix2=0,nx_kr-1	!x
 do iy2=0,ny_kr-1	!y
  do iz2=0,nz-1	!z
   do icf2=0,1	!c/f
    do is2=0,1	!up/down 
     do o2=0,l_orb(icf2)-1
      ind2p=index_xyzcso (ix2 , iy2 , iz2, icf2, is2, o2, nx_kr, ny_kr, nz)

!non diag ke (H_0 and perturbation)
if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*(-disp_n+ky2)	!td12
if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*(-disp_n+ky2)	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*(-disp_n+ky2)	!td12
if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*(-disp_n+ky2)	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78
if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78



!kinetic energy for discretized model --> iz is layer index
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tc_*g1d	!d1
!if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tc_*g2d	!d2
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g1f	!f1
!if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g2f	!f2
!if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==2 .and. is==is2 .and. (iz2==iz+1 .or. iz2==iz-1)) ham(indp,ind2p)=tf_*b_k**2*g7	!f7

!hybridization for discretized model --> iz is layer index
!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=+f1*(2*is-1)*Vcf_*b_k/2!*pi		!df1
!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=-f1*(2*is-1)*Vcf_*b_k/2!*pi		!df1
!if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=+f2*(2*is-1)*Vcf_*b_k/2!*pi		!df2
!if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=-f2*(2*is-1)*Vcf_*b_k/2!*pi		!df2
!if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=-f1*(2*is-1)*Vcf_*b_k/2!*pi		!fd1
!if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=+f1*(2*is-1)*Vcf_*b_k/2!*pi		!fd1
!if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz+1) ham(indp,ind2p)=-f2*(2*is-1)*Vcf_*b_k/2!*pi		!fd2
!if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz-1) ham(indp,ind2p)=+f2*(2*is-1)*Vcf_*b_k/2!*pi		!fd2



!if (iz.ne.iz2)sincos= 2*(-1+(-1)**(iz+iz2))*(iz2+1)*(iz+1)/((iz+1)**2-(iz2+1)**2)/(nz+1)
!if (iz.ne.iz2)sincos= 2*(-1+(-1)**(iz+iz2))*sin(pi*(iz2+1)*(iz+1)/(nz+1))/((iz+1)**2-(iz2+1)**2)/pi
!if (iz.ne.iz2)sincos=2*sin(pi*(iz+1)*(iz2+1)/(nz+1))	!discretized hybr (try)
if (iz.ne.iz2) sincos=hybr_k(iz2+1,iz+1,nz)*norm_k(iz+1,nz)*norm_k(iz2+1,nz)

!hybr kx-> -i dx (belongs to H_0)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-h1*sincos*Vcf_*b_k!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-h2*sincos*Vcf_*b_k!*2*pi		!df2
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-hx*sincos*Vcf_*b_k			!df12
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-hx*sincos*Vcf_*b_k			!df21
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=h1*sincos*Vcf_*b_k!*2*pi		!fd1
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=h2*sincos*Vcf_*b_k!*2*pi		!fd2
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=hx*sincos*Vcf_*b_k			!fd12
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=hx*sincos*Vcf_*b_k			!fd21
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-h7*sincos*Vcf_*b_k!*2*pi		!df17
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=-h72*sincos*Vcf_*b_k			!df27
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==0 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=h7*sincos*Vcf_*b_k!*2*pi		!fd71
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==1 .and. is.ne.is2 .and. iz2.ne.iz) ham(indp,ind2p)=h72*sincos*Vcf_*b_k			!fd72


!if (iz.ne.iz2)sincos2=2*(-1+(-1)**(iz2+iz))*(iz+1)*(iz2+1)/((iz2+1)**2-(iz+1)**2)/(nz+1)
!if (iz.ne.iz2)sincos2=-sincos

!kin energy t12 k^2_y-k^2_x, no dz term
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf2==0 .and. icf==0 .and. o2==0 .and. o1==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf2==0 .and. icf==0 .and. o2==1 .and. o1==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*k2xy
!if (icf2==1 .and. icf==1 .and. o2==0 .and. o1==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!if (icf2==1 .and. icf==1 .and. o2==1 .and. o1==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*k2xy*b_k**2
!hybr ky (perturbation)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h1*km	!df
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h1*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*kp
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*hx*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h2*km
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h2*kp
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h7*km	!df7
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h7*kp
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h72*kp	!df72
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*h72*km

if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h1*kp	!fd
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h1*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*km
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*hx*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h2*kp
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h2*km
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h7*kp	!fd7
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h7*km
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h72*km	!fd72
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=(0,1)*h72*kp
!hybr kz(perturbation)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*f1*kz*(2*is-1)	!df
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*f2*kz*(2*is-1)	!df
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=-(0,1)*f7*kz*(2*is-1)	!df
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=(0,1)*f1*kz*(2*is-1)	!fd
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=(0,1)*f2*kz*(2*is-1)	!fd
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=(0,1)*f7*kz*(2*is-1)	!fd

     enddo
    enddo
   enddo
  enddo
 enddo
enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo




do i=0, size_h_kr-1
do j=i, size_h_kr-1
!if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9 .and. abs(Ham(j,i)+conjg(Ham(i,j)))>1e-9) then
if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham(i,j), Ham(j,i)
!stop
endif
enddo
enddo

end subroutine build_ham_kr_kdotp_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_ham_kr_kdotp_theta(ik,Ham,b, lambda,mu,k_vecs,nk,ldiag,theta,phi,omega,kx0,ky0)
complex(KIND=idpc),intent(out) 	:: Ham(0:size_h_kr-1,0:size_h_kr-1)
real(KIND=idp),intent(in) 	:: b(0:n_sites_kr-1), lambda(0:n_sites_kr-1), mu,k_vecs(0:nk-1,6)
integer, intent(in)::ik,nk
logical, intent(in):: ldiag
real(KIND=idp)::disp_n,ky2,kz2,kz,theta,phi,omega,k_kbar(3,3),kbar_x,kbar_y
real(KIND=idp)::k2x_xx,k2x_yy,k2x_zz,k2x_xy,k2x_xz,k2x_yz,k2y_xx,k2y_yy,k2y_zz,k2y_xy,k2y_xz,k2y_yz,k2z_xx,k2z_yy,k2z_zz,k2z_xy,k2z_xz,k2z_yz,k_xx,k_yy,k_zxzx,k_zxxz
real(KIND=idp)::k2parm_par, k2parm_z,kx0,ky0,k2par_par,k2z_par,k2par_z,e_test(0:2*nz-1)
complex(KIND=idpc)::kp, km,kp_x,kp_y,kp_z, km_x,km_y,km_z,kz_x,kz_y,kz_z,sincos,sincos2, test_m(0:2*nz-1,0:2*nz-1), w_test(0:2*nz-1,0:2*nz-1)

!bar k_z -> -i d/dz

!print *, "build ham kr kdp"
ham=0

!k as a fct of kbar -- rotk1 in mathematica
k_kbar(1,1)=cos(omega)*cos(theta)*cos(phi)-sin(omega)*sin(phi)
k_kbar(1,2)=-cos(omega)*sin(phi)-sin(omega)*cos(theta)*cos(phi)
k_kbar(1,3)=sin(theta)*cos(phi)
k_kbar(2,1)=cos(omega)*cos(theta)*sin(phi)+sin(omega)*cos(phi)
k_kbar(2,2)=cos(omega)*cos(phi)-sin(omega)*cos(theta)*sin(phi)
k_kbar(2,3)=sin(theta)*sin(phi)
k_kbar(3,1)=-cos(omega)*sin(theta)
k_kbar(3,2)=sin(omega)*sin(theta)
k_kbar(3,3)=cos(theta)

!print *, "rot mat k(kbar)"
!do i=1,3
!write (*, "(3f10.5)") (k_kbar(i,j),j=1,3)
!enddo

!kbar close to 0

kbar_x=(k_vecs(ik,1))
kbar_y=(k_vecs(ik,2))

if (surf_type=="110") kbar_y=kbar_y/sqrt(2.) 
if (surf_type=="210") kbar_y=kbar_y/sqrt(5.) 

kbar_x=kbar_x-kx0
kbar_y=kbar_y-ky0

write (*, "(10f10.5)"), k_vecs(ik,1), k_vecs(ik,2), kbar_x, kbar_y,kx0,ky0



kbar_x=kbar_x*2*pi
kbar_y=kbar_y*2*pi


!kp=+(0,1)*k_vecs(ik,2)*Vcf_*b_k*2*pi	!+iky
!km=-(0,1)*k_vecs(ik,2)*Vcf_*b_k*2*pi	!-iky
!kz=(k_vecs(ik,1)-0.5)*Vcf_*b_k*2*pi	!kz
!ky2=k_vecs(ik,2)**2*4*pi**2		!ky^2
!kz2=(k_vecs(ik,1)-0.5)**2*4*pi**2	!kz^2

!k^2 as a fct of kbar^2
k2x_xx=k_kbar(1,1)**2*kbar_x**2
k2x_yy=k_kbar(1,2)**2*kbar_y**2
k2x_zz=k_kbar(1,3)**2
k2x_xy=k_kbar(1,1)*k_kbar(1,2)*kbar_y*kbar_x
k2x_yz=k_kbar(1,2)*k_kbar(1,3)*kbar_y
k2x_xz=k_kbar(1,3)*k_kbar(1,1)*kbar_x
k2y_xx=k_kbar(2,1)**2*kbar_x**2
k2y_yy=k_kbar(2,2)**2*kbar_y**2
k2y_zz=k_kbar(2,3)**2
k2y_xy=k_kbar(2,1)*k_kbar(2,2)*kbar_y*kbar_x
k2y_yz=k_kbar(2,2)*k_kbar(2,3)*kbar_y
k2y_xz=k_kbar(2,3)*k_kbar(2,1)*kbar_x
k2z_xx=k_kbar(3,1)**2*kbar_x**2
k2z_yy=k_kbar(3,2)**2*kbar_y**2
k2z_zz=k_kbar(3,3)**2
k2z_xy=k_kbar(3,1)*k_kbar(3,2)*kbar_y*kbar_x
k2z_yz=k_kbar(3,2)*k_kbar(3,3)*kbar_y
k2z_xz=k_kbar(3,3)*k_kbar(3,1)*kbar_x


k2par_par=k2x_xx+k2x_yy+k2x_xy+k2y_xx+k2y_yy+k2y_xy
k2z_par=k2z_xx+k2z_yy+k2z_xy
k2par_z=k2x_zz+k2y_zz

!kzyzy=k_kbar(3,3)**2*k_kbar(2,2)
k_zxzx=k_kbar(3,3)**2 !*k_kbar(1,1)	!only for theta=0
k_zxxz=k_kbar(1,3)**2 !*k_kbar(3,1)	!only for theta=pi/2


!here iz denotes n, the index of the wfc sin [pi n z/(N+1)], n,z=1,..,N


do ix=0,nx_kr-1	!x
 do iy=0,ny_kr-1	!y
  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1
          
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx_kr, ny_kr, nz)

!y^disp_n=(iz+1.)**2/(nz+1.)**2*pi**2	!free dispersion
if (.not. pbcz) disp_n=2-2*cos(pi*(iz+1)/(nz+1))	!discretized dispersion
if (      pbcz) disp_n=2-2*cos(2.*pi*iz/nz)	!discretized dispersion

!onsite enrgy, kx^2+ky^2, dz^2
!iz k index of sin(pi*x*k)/(N+1)
!									H_0				|	H_P (k^2)
if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*disp_n*(l1d*k2par_z+g1d*k2z_zz)			+tc_*(l1d*k2par_par+g1d*k2z_par)
if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*disp_n*(l2d*k2par_z+g2d*k2z_zz)			+tc_*(l2d*k2par_par+g2d*k2z_par)
if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*disp_n*(l1f*k2par_z+g1f*k2z_zz)	+tf_*b_k**2*(l1f*k2par_par+g1f*k2z_par)
if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*disp_n*(l2f*k2par_z+g2f*k2z_zz)	+tf_*b_k**2*(l2f*k2par_par+g2f*k2z_par)
if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tf_*b_k**2*disp_n*(l7 *k2par_z+g7 *k2z_zz)	+tf_*b_k**2*(l7 *k2par_par+g7 *k2z_par)

!if (icf==0 .and. o1==0) test_m(2*iz,2*iz)=eckp1-mu_k+tc_*disp_n*(l1d*k2par_z+g1d*k2z_zz)
!if (icf==1 .and. o1==0) test_m(2*iz+1,2*iz+1)=efkp1-mu_k-lambda_k+tf_*b_k**2*disp_n*(l1f*k2par_z+g1f*k2z_zz)

!if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tf_*b_k**2*kpar2*l7 +tf_*b_k**2*g7*disp_n
!iz layer index
!if (icf==0 .and. o1==0) ham(indp,indp)=eckp1-mu_k+tc_*kpar2*l1d+2*tc_*g1d
!if (icf==0 .and. o1==1) ham(indp,indp)=eckp2-mu_k+tc_*kpar2*l2d+2*tc_*g2d
!if (icf==1 .and. o1==0) ham(indp,indp)=efkp1-mu_k-lambda_k+tf_*b_k**2*kpar2*l1f+2*tf_*b_k**2*g1f
!if (icf==1 .and. o1==1) ham(indp,indp)=efkp2-mu_k-lambda_k+tf_*b_k**2*kpar2*l2f+2*tf_*b_k**2*g2f
!if (icf==1 .and. o1==2) ham(indp,indp)=efkp7-mu_k-lambda_k+tc_*b_k**2*kpar2*l7+ 2*tf_*b_k**2*g7
!m78


do iz2=0,nz-1	!z
!if (iz.ne.iz2)sincos= 2*(-1.+(-1.)**(iz+iz2))*(iz2+1.)*(iz+1.)/((iz+1.)**2-(iz2+1.)**2)/(nz+1.)	!continous
if (.not. pbcz) sincos=hybr_k(iz2+1,iz+1,nz)*norm_k(iz+1,nz)*norm_k(iz2+1,nz)			!discretized
if (      pbcz) sincos=hybr_k_pbc(iz2,iz,nz)
!if (abs(sincos)>1e-5) print *, iz,iz2,sincos
!sincos2=-sincos

 do ix2=0,nx_kr-1	!x
  do iy2=0,ny_kr-1	!y
   do icf2=0,1	!c/f
    do is2=0,1	!up/down 
     do o2=0,l_orb(icf2)-1
      ind2p=index_xyzcso (ix2 , iy2 , iz2, icf2, is2, o2, nx_kr, ny_kr, nz)

!non diag ke (H_0 and perturbation)
k2parm_par=-k2x_xx-k2x_yy-k2x_xy+k2y_xx+k2y_yy+k2y_xy
k2parm_z=-k2x_zz+k2y_zz

if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*(disp_n*k2parm_z+k2parm_par)	!td12
if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tc_*td12*(disp_n*k2parm_z+k2parm_par)	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*(disp_n*k2parm_z+k2parm_par)	!td12
if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=tf_*tf12*(disp_n*k2parm_z+k2parm_par)	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78
if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=m78	!m78


!kinetic energy terms in bar kx * bar kz and bar ky * bar kz
!diagonal td tf
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tc_*(l1d*(k2x_yz+k2x_xz+k2y_yz+k2y_xz)+g1d*(k2z_xz+k2z_yz))*sincos
!if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tc_*(l2d*(k2x_yz+k2x_xz+k2y_yz+k2y_xz)+g2d*(k2z_xz+k2z_yz))*sincos
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tf_*b_k**2*(l1f*(k2x_yz+k2x_xz+k2y_yz+k2y_xz)+g1f*(k2z_xz+k2z_yz))*sincos
!if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tf_*b_k**2*(l2f*(k2x_yz+k2x_xz+k2y_yz+k2y_xz)+g2f*(k2z_xz+k2z_yz))*sincos
!if (icf==1 .and. icf2==1 .and. o1==2 .and. o2==2 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tf_*b_k**2*(l7 *(k2x_yz+k2x_xz+k2y_yz+k2y_xz)+g7 *(k2z_xz+k2z_yz))*sincos

!nondiagonal td12 tf12
if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tc_*td12*sincos*(k2y_yz+k2y_xz-k2x_yz+k2x_xz)	!td12
if (icf==0 .and. icf2==0 .and. o1==1 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tc_*td12*sincos*(k2y_yz+k2y_xz-k2x_yz+k2x_xz)	!td12
if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tf_*b_k**2*tf12*sincos*(k2y_yz+k2y_xz-k2x_yz+k2x_xz)	!tf12
if (icf==1 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*tf_*b_k**2*tf12*sincos*(k2y_yz+k2y_xz-k2x_yz+k2x_xz)	!tf12


kp_x=(k_kbar(1,1)+(0,1)*k_kbar(2,1))*kbar_x
kp_y=(k_kbar(1,2)+(0,1)*k_kbar(2,2))*kbar_y
kp_z=(k_kbar(1,3)+(0,1)*k_kbar(2,3))
km_x=(k_kbar(1,1)-(0,1)*k_kbar(2,1))*kbar_x
km_y=(k_kbar(1,2)-(0,1)*k_kbar(2,2))*kbar_y
km_z=(k_kbar(1,3)-(0,1)*k_kbar(2,3))
kz_x=k_kbar(3,1)*kbar_x
kz_y=k_kbar(3,2)*kbar_y
kz_z=k_kbar(3,3)

kp_x=(k_kbar(1,1)+(0,1)*k_kbar(2,1))*kbar_x


k_xx=k_kbar(1,1)*kbar_x
k_yy=k_kbar(2,2)*kbar_y

!print *, "kz_z=",kz_z

!hybr kp km bar kz-> -i dz (belongs to H_0)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-h1*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-h1*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-h2*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-h2*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-hx*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-hx*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-hx*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-hx*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-h7*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-h7*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==1 .and. is2==0) ham(indp,ind2p)=ham(indp,ind2p)-h72*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==0 .and. is2==1) ham(indp,ind2p)=ham(indp,ind2p)-h72*sincos*Vcf_*b_k*km_z!*2*pi		!df1

if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==0 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+h1*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==0 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+h1*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==1 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+h2*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==1 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+h2*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==1 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+hx*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==1 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+hx*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==0 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+hx*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==0 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+hx*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==2 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+h7*sincos*Vcf_*b_k*kp_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==0 .and. o1==2 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+h7*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==2 .and. is2==1 .and. is==0) ham(indp,ind2p)=ham(indp,ind2p)+h72*sincos*Vcf_*b_k*km_z!*2*pi		!df1
if (icf==1 .and. icf2==0 .and. o2==1 .and. o1==2 .and. is2==0 .and. is==1) ham(indp,ind2p)=ham(indp,ind2p)+h72*sincos*Vcf_*b_k*kp_z!*2*pi		!df1

!hybr kz bar kz-> -i dz (belongs to H_0)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-f1*sincos*(2*is-1)*Vcf_*b_k*kz_z		!df1
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-f2*sincos*(2*is-1)*Vcf_*b_k*kz_z		!df2
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)-f7*sincos*(2*is-1)*Vcf_*b_k*kz_z		!df7
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)+f1*sincos*(2*is-1)*Vcf_*b_k*kz_z		!fd1
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)+f2*sincos*(2*is-1)*Vcf_*b_k*kz_z		!fd2
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==0 .and. is==is2) ham(indp,ind2p)=ham(indp,ind2p)+f7*sincos*(2*is-1)*Vcf_*b_k*kz_z		!fd7


!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) print *, ham(indp,ind2p)
!if (icf==0 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) print *, ham(indp,ind2p)
!if (icf==1 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) print *, ham(indp,ind2p)

!if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz==iz2) test_m(2*iz,2*iz+1)=-f1*sincos*Vcf_*b_k*kz_z!*2*pi		!df1
!if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz==iz2) test_m(2*iz+1,2*iz)=+f1*sincos*Vcf_*b_k*kz_z!*2*pi		!df1



!hybr kx ky kbar x kbar y (perturbation)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*(h1+h1k2*disp_n*k_zxzx)*Vcf_*b_k*(km_x+km_y)	!df
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*(h1+h1k2*disp_n*k_zxzx)*Vcf_*b_k*(kp_x+kp_y)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*hx*Vcf_*b_k*(kp_x+kp_y)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*hx*Vcf_*b_k*(km_x+km_y)
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*hx*Vcf_*b_k*(kp_x+kp_y)
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==0 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*hx*Vcf_*b_k*(km_x+km_y)
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h2*Vcf_*b_k*(km_x+km_y)
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h2*Vcf_*b_k*(kp_x+kp_y)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h7*Vcf_*b_k*(km_x+km_y)	!df7
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h7*Vcf_*b_k*(kp_x+kp_y)
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==1 .and. is2==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h72*Vcf_*b_k*(kp_x+kp_y)	!df72
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==2 .and. is==0 .and. is2==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*h72*Vcf_*b_k*(km_x+km_y)

if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*(h1+h1k2*disp_n*k_zxzx)*Vcf_*b_k*(kp_x+kp_y)	!fd
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*(h1+h1k2*disp_n*k_zxzx)*Vcf_*b_k*(km_x+km_y)
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*hx*Vcf_*b_k*(km_x+km_y)
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*hx*Vcf_*b_k*(kp_x+kp_y)
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*hx*Vcf_*b_k*(km_x+km_y)
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==0 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*hx*Vcf_*b_k*(kp_x+kp_y)
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h2*Vcf_*b_k*(kp_x+kp_y)
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==1 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h2*Vcf_*b_k*(km_x+km_y)
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h7*Vcf_*b_k*(kp_x+kp_y)	!fd7
if (icf2==0 .and. icf==1 .and. o2==0 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h7*Vcf_*b_k*(km_x+km_y)
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==1 .and. is==0 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h72*Vcf_*b_k*(km_x+km_y)	!fd72
if (icf2==0 .and. icf==1 .and. o2==1 .and. o1==2 .and. is2==0 .and. is==1 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*h72*Vcf_*b_k*(kp_x+kp_y)
!hybr kz kbar x kbar y(perturbation)
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*f1*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!df
if (icf==0 .and. icf2==1 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*f2*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!df
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==2 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*f7*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!df
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*f1*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!fd
if (icf==1 .and. icf2==0 .and. o1==1 .and. o2==1 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*f2*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!fd
if (icf==1 .and. icf2==0 .and. o1==2 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*f7*(2*is-1)*Vcf_*b_k*(kz_x+kz_y)	!fd
!correction from 3rd order terms for theta=pi/2
if (icf==0 .and. icf2==1 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)-(0,1)*f1k2*(2*is-1)*Vcf_*b_k*kz_x*disp_n*k_zxxz	!df
if (icf==1 .and. icf2==0 .and. o1==0 .and. o2==0 .and. is==is2 .and. iz2==iz) ham(indp,ind2p)=ham(indp,ind2p)+(0,1)*f1k2*(2*is-1)*Vcf_*b_k*kz_x*disp_n*k_zxxz	!df

     enddo
    enddo
   enddo
  enddo
 enddo
enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo

!print *, "test_m"

!CALL DIAGONALIZE( 'V', 2*nz, 2*nz, test_m, E_test, W_test )

!do i=0,2*nz-1
!print *, e_test(i)
!enddo

!if (ik==141) then
!OPEN(unit=231,file=trim(label)//'/test_m',status='unknown')

!do i=0,2*nz-1
!write(231, "(i5,f15.7)"), ik,e_test(i)
!enddo

! close(unit=231)
!endif


do i=0, size_h_kr-1
do j=i, size_h_kr-1
!if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9 .and. abs(Ham(j,i)+conjg(Ham(i,j)))>1e-9) then
if (abs(Ham(j,i)-conjg(Ham(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham(i,j), Ham(j,i)
!stop
endif
enddo
enddo

end subroutine build_ham_kr_kdotp_theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine build_ham_kr_read(ik,Ham_kr,b_kr,lambda_kr, mu_kr, ldiag,nk,k_vecs)
complex(KIND=idpc),intent(out) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1)
real(KIND=idp),intent(in) 	:: b_kr(0:n_sites_kr-1), lambda_kr(0:n_sites_kr-1), mu_kr,k_vecs(0:nk-1,6)
integer, intent(in)		::ik,nk
integer				:: i
logical, intent(in)		::ldiag


Ham_kr=0


do iz=0, nz-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
    do iz2=0, nz-1
     do icf2=0,1
      do is2=0,1
       do o2=0,l_orb(icf2)-1
 
 indp=index_zcso (iz, icf, is, o1, nz)
 ind2p=index_zcso (iz2, icf2, is2, o2, nz)
 i=index_cso (icf, is, o1)
 j=index_cso (icf2, is2, o2)


do lx=-n_shells,n_shells
do ly=-n_shells,n_shells
do lz=-n_shells,n_shells

if (sqrt(real(lx**2+ly**2+lz**2)) <=(dist_max+0.001) .and. sqrt(real(lx**2+ly**2+lz**2)) >0.001 ) then !I exclude on site and too far shells

 if((iz==iz2+lz .and. .not. pbcz) .or. (mod(iz-iz2-lz+6*nz,nz)==0 .and. pbcz)) then

!dd
 if (icf==0 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly)*(0,1))

!df

if (icf==1 .and. icf2==0) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly)*(0,1))*b_kr(iz)
if (icf==0 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly)*(0,1))*b_kr(iz2)

!ff
if (icf==1 .and. icf2==1) Ham_kr(indp, ind2p)=Ham_kr(indp, ind2p)+Ham_input(lx,ly,lz,i,j)*exp(2*pi*(k_vecs(ik,1)*lx+k_vecs(ik,2)*ly)*(0,1))*b_kr(iz)*b_kr(iz2)


 endif
endif
enddo
enddo
enddo

!on site energy --> I only consider energy, not non-diag elements, which should be zero
if (ldiag .and. indp==ind2p) then

 if (icf==0) Ham_kr(indp,indp)=Ham_kr(indp,indp)+Ham_input(0,0,0,i,i)-mu_kr+ec_site_kr(iz,o1)			!lambda and mu
 if (icf==1) Ham_kr(indp,indp)=Ham_kr(indp,indp)+Ham_input(0,0,0,i,i)-mu_kr-lambda_kr(iz)+ef_site_kr(iz,o1)

endif

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo




!if (n_sites <5) then
do i=0, size_h_kr-1
do j=i, size_h_kr-1
if (abs(Ham_kr(j,i)-conjg(Ham_kr(i,j)))>1e-9) then
print *, "error in Hamiltonian kr", i, j, Ham_kr(i,j), Ham_kr(j,i)
stop
endif
enddo
enddo
!endif

end subroutine build_ham_kr_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_ham_kr_2d(ik,Ham_kr,nk,k_vecs)
complex(KIND=idpc),intent(out) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1)
integer, intent(in)		::ik,nk
integer				:: i
real(KIND=idp),intent(in) 	:: k_vecs(0:nk-1,6)


Ham_kr=0



!zones with circles
!if (((k_vecs(ik,1)-0.5)**2+k_vecs(ik,2)**2)<1./8.) then
!Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)-0.5)*v_1	!pi,0
!elseif (((k_vecs(ik,1)+0.5)**2+k_vecs(ik,2)**2)<1./8.) then
!Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)+0.5)*v_1	!-pi,0
!elseif ((k_vecs(ik,1)**2+(k_vecs(ik,2)-0.5)**2)<1./8.) then
!Ham_kr(0,1)=(k_vecs(ik,2)-0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,pi
!elseif ((k_vecs(ik,1)**2+(k_vecs(ik,2)+0.5)**2)<1./8.) then
!Ham_kr(0,1)=(k_vecs(ik,2)+0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,-pi
!else
!Ham_kr(0,1)=1e5
!endif

!zones with lines
if (k_vecs(ik,1)>=0 .and. abs(k_vecs(ik,1))>=abs(k_vecs(ik,2))) then
Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)-0.5)*v_1	!pi,0
elseif (k_vecs(ik,1)<=0 .and. abs(k_vecs(ik,1))>=abs(k_vecs(ik,2))) then
Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)+0.5)*v_1	!-pi,0
elseif (k_vecs(ik,2)>=0 .and. abs(k_vecs(ik,2))>=abs(k_vecs(ik,1))) then
Ham_kr(0,1)=(k_vecs(ik,2)-0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,pi
elseif (k_vecs(ik,2)<=0 .and. abs(k_vecs(ik,2))>=abs(k_vecs(ik,1))) then
Ham_kr(0,1)=(k_vecs(ik,2)+0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,-pi
else
Ham_kr(0,1)=1e5
endif

!energy cutoff
!if (k_vecs(ik,2)**2*v_2**2+(k_vecs(ik,1)-0.5)**2*v_1**2<omega_cutoff**2) then
!Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)-0.5)*v_1	!pi,0
!elseif (k_vecs(ik,2)**2*v_2**2+(k_vecs(ik,1)+0.5)**2*v_1**2<omega_cutoff**2) then
!Ham_kr(0,1)=k_vecs(ik,2)*v_2-(0,1)*(k_vecs(ik,1)+0.5)*v_1	!-pi,0
!elseif ((k_vecs(ik,2)-0.5)**2*v_1**2+k_vecs(ik,1)**2*v_2**2<omega_cutoff**2) then
!Ham_kr(0,1)=(k_vecs(ik,2)-0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,pi
!elseif ((k_vecs(ik,2)+0.5)**2*v_1**2+k_vecs(ik,1)**2*v_2**2<omega_cutoff**2) then
!Ham_kr(0,1)=(k_vecs(ik,2)+0.5)*v_1-(0,1)*k_vecs(ik,1)*v_2	!0,-pi
!else
!Ham_kr(0,1)=1e5
!endif



Ham_kr(1,0)=conjg(Ham_kr(0,1))

end subroutine build_ham_kr_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagonalize_ham_kr(Ham_kr, E_kr, W_kr)
#ifdef cmplx
complex(KIND=idpc) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1),W_kr(0:size_h_kr-1,0:size_h_kr-1)
#else
real(KIND=idpc) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1),W_kr(0:size_h_kr-1,0:size_h_kr-1)
#endif
real(KIND=idp)		:: E_kr(0:size_h_kr-1)

CALL DIAGONALIZE( 'V', size_h_kr, size_h_kr, Ham_kr, E_kr, W_kr )


end subroutine diagonalize_ham_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagonalize_ham_kr_kdotp(Ham_kr, E_kr, W_kr)
#ifdef cmplx
complex(KIND=idpc) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1),W_kr(0:size_h_kr-1,0:size_h_kr-1)
#else
real(KIND=idpc) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1),W_kr(0:size_h_kr-1,0:size_h_kr-1)
#endif
complex(KIND=idpc)	:: E_kr(0:size_h_kr-1)

CALL DIAGONALIZE_s( size_h_kr, Ham_kr, E_kr)


end subroutine diagonalize_ham_kr_kdotp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagonalize_ham_kr2(Ham_kr, E_kr, W_kr)
#ifdef cmplx
complex(KIND=idpc) 	:: Ham_kr(0:size_h_kr/2-1,0:size_h_kr/2-1),W_kr(0:size_h_kr/2-1,0:size_h_kr/2-1)
#else
real(KIND=idpc) 	:: Ham_kr(0:size_h_kr/2-1,0:size_h_kr/2-1),W_kr(0:size_h_kr/2-1,0:size_h_kr/2-1)
#endif
real(KIND=idp)		:: E_kr(0:size_h_kr/2-1)

CALL DIAGONALIZE( 'V', size_h_kr/2, size_h_kr/2, Ham_kr, E_kr, W_kr )


end subroutine diagonalize_ham_kr2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_mu_kr (e_kr_tot)
real(KIND=idp), intent(in)	:: E_kr_tot(0:size_h_kr-1,0:nk2d-1)
integer			        :: i,j,ik
real (kind=idp)			::thresh, alfa_mu, mu_old_kr


mu_old_kr=mu_kr

i_iter_mu=0
diff_mu=100
alfa_mu=alpha_mu_kr

write(*,"(a40,100f15.10)") "******* Mu loop, starting mu=", mu_kr
!write(*,"(a5,5(a20),a20,a15)") "i_mu","mu",  "Nfctot",  "N_el", "diff_mu", "mu_new", "dndmu"

!!!!!!!mu loop


!	if (i_iter<5) then
!	thresh=0.1
!	alfa_mu=alpha_mu_kr*0.001
!	else
!	thresh=thr_mu
	alfa_mu=alpha_mu_kr
!	endif

!alfa_mu=0.01*(1-1/(i_iter+0.1))
!thresh=thr_mu+0.01/i_iter

!alfa_mu=alpha_mu_k
thresh=thr_mu

mu_old_kr=mu_kr
mu_new_kr=mu_kr

	do while (diff_mu>thresh)
	i_iter_mu=i_iter_mu+1

nfctot=0
   
   do i=0, size_h_kr-1	!eigenvalues
    do ik=0,nk2d-1	!k points
!    nfctot=nfctot+1/(exp((e_kr_tot(i,ik)-mu_kr)/T)+1)*k_vecs2d(ik,4)
    nfctot=nfctot+1/(exp((e_kr_tot(i,ik)+mu_old_kr-mu_new_kr)/T)+1)*k_vecs2d(ik,4)
   enddo
   enddo
!!!!!!Difference between wanted and obtained N_el

diff_mu=abs(N_el_kr-Nfctot)

!!!!!!!!!!!!! Change mu

dndmu=0
do i=0,size_h_kr-1
do ik=0, nk2d-1
dndmu=dndmu+1/T/(exp((e_kr_tot(i,ik)+mu_old_kr-mu_new_kr)/T)+1)/(exp((-e_kr_tot(i,ik)+mu_old_kr-mu_new_kr)/T)+1)*k_vecs2d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
enddo
enddo




!if (abs(N_el_k-Nfctot)<N_k*0.1)then
dndmu=max(dndmu,0.001)		!to avoid dndmu --> 0 when T-->0
!else
!dndmu=-N_k/(4*tc_)
!endif

mu_new_kr=mu_kr+(N_el_kr-Nfctot)/dndmu

!print *, dndmu, mu_new_k

	if (mod(i_iter_mu,10001)==0) write(*,"(i6,5f20.10,f20.10, f15.10)") i_iter_mu, mu_kr,  Nfctot,  N_el_kr, diff_mu, mu_new_kr, dndmu


	if (diff_mu>0.5) then
!	thresh=0.1
	alfa_mu=alpha_mu_kr*0.01
	else
!	thresh=thr_mu
	alfa_mu=alpha_mu_kr
	endif


	mu_new_kr=alfa_mu*mu_new_kr+(1-alfa_mu)*mu_kr
	mu_kr=mu_new_kr

	enddo


!!!!!!!!!!Finish loop



write(*,"(a60,f15.10, i10, f15.5)") "********Mu loop finished, final mu, n_iter, n_fctot=", mu_new_kr, i_iter_mu, Nfctot

!pause
end subroutine compute_mu_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_mu_kr_1 (e_kr_tot)
real(KIND=idp), intent(in)	:: E_kr_tot(0:size_h_kr-1,0:nk2d-1)
integer			        :: i,j,ik
real (kind=idp)			::thresh, alfa_mu, mu_old_kr


!mu_old_kr=mu_kr

!i_iter_mu=0
!diff_mu=100
!alfa_mu=alpha_mu_kr

!write(*,"(a40,100f15.10)") "******* Mu loop, starting mu=", mu_kr
!write(*,"(a5,5(a20),a20,a15)") "i_mu","mu",  "Nfctot",  "N_el", "diff_mu", "mu_new", "dndmu"

!!!!!!!mu loop


!	if (i_iter<5) then
!	thresh=0.1
!	alfa_mu=alpha_mu_kr*0.001
!	else
!	thresh=thr_mu
!	alfa_mu=alpha_mu_kr
!	endif

!alfa_mu=0.01*(1-1/(i_iter+0.1))
!thresh=thr_mu+0.01/i_iter

!alfa_mu=alpha_mu_k
!thresh=thr_mu



!	do while (diff_mu>thresh)
!	i_iter_mu=i_iter_mu+1

nfctot=0
   
   do i=0, size_h_kr-1	!eigenvalues
    do ik=0,nk2d-1	!k points
    nfctot=nfctot+1/(exp((e_kr_tot(i,ik)-mu_kr)/T)+1)*k_vecs2d(ik,4)
   enddo
   enddo
!!!!!!Difference between wanted and obtained N_el

diff_mu=abs(N_el_kr-Nfctot)

!!!!!!!!!!!!! Change mu

dndmu=0
do i=0,size_h_kr-1
do ik=0, nk2d-1
dndmu=dndmu+1/T/(exp((e_kr_tot(i,ik)-mu_kr)/T)+1)/(exp((-e_kr_tot(i,ik)+mu_kr)/T)+1)*k_vecs2d(ik,4)		!d(n)/d(mu)=sum_i f'(e_i-mu)
enddo
enddo




!if (abs(N_el_k-Nfctot)<N_k*0.1)then
dndmu=max(dndmu,0.001)		!to avoid dndmu --> 0 when T-->0
!else
!dndmu=-N_k/(4*tc_)
!endif

mu_new_kr=mu_kr+(N_el_kr-Nfctot)/dndmu
mu_new_kr=alpha_mu_kr*mu_new_kr+(1-alpha_mu_kr)*mu_kr

!print *, dndmu, mu_new_k

!	if (mod(i_iter_mu,10001)==0) write(*,"(i6,5f20.10,f20.10, f15.10)") i_iter_mu, mu_kr,  Nfctot,  N_el_kr, diff_mu, mu_new_kr, dndmu


!	if (diff_mu>0.5) then
!	thresh=0.1
!	alfa_mu=alpha_mu_kr*0.0001
!	else
!thresh=thr_mu
!	alfa_mu=alpha_mu_kr
!	endif


!	mu_new_kr=alfa_mu*mu_new_kr+(1-alfa_mu)*mu_kr
!	mu_kr=mu_new_kr

!	enddo


!!!!!!!!!!Finish loop

!mu_kr=mu_old_kr


write(*,"(a60,2f15.10, i10, f15.5)") "Mu old mu,final mu, n_iter, n_fctot=", mu_kr,mu_new_kr, i_iter_mu, Nfctot

!pause
end subroutine compute_mu_kr_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_nf_kr (nf_kr,nc_kr,b_new_kr,e_kr_tot,w_kr_tot,ef_kr, lambda_new_kr, nfc_kr,diff,nff_kr,mu)
complex(KIND=idpc) 	:: w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:n_k_kr-1)
real(KIND=idp) 		:: b_new_kr(0:n_sites_kr-1),lambda_new_kr(0:n_sites_kr-1)
real(KIND=idp) 		:: E_kr_tot(0:size_h_kr-1,0:n_k_kr-1),diff,mu
real(kind=idp)		:: nf_kr(0:n_sites_kr-1),nc_kr(0:n_sites_kr-1),nfc_kr(0:n_sites_kr-1),nff_kr(0:n_sites_kr-1)
real(kind=idp)		:: ec_kr(0:n_sites_kr-1),ef_kr(0:n_sites_kr-1)
integer		::i,j,ik,is,is2

!only works for nx_kr=1 ny_kr=1

nf_kr=0		!sum_k <f^+f>
nfk_kr=0	!sum_k tf*f_kr*<f^+f>
nc_kr=0		!sum_k <c^+c>
nfck_kr=0	!sum_k Phi_k <c_kz^+f_kz>
nfc_kr=0	!sum_kz' V Phi_zz' <c_kz^+f_kz'>
nff_kr=0	!sum_kz' tf*b_z'* <f_kz^+f_kz'>
  
  
do ik=0,nk2d-1		!k -points
!print *, ik
! call factor_kr(ik,f_kr)
!  f_kr=factor_kr_tot(ik)
! call build_ham_kr(ik,V_kr,one_kr,zero_kr, 0.0d0,.false.,factor_kr_tot,nk2d,phi_kr_sigma_alpha_tot,phi_kr_alpha_sigma_tot)		!da fare v_kr_tot


 do i=0,size_h_kr-1		!eigenvalues

  occup=1/(exp((e_kr_tot(i,ik)-mu)/T)+1) * k_vecs2d(ik,4)		!mu=0 already in hamiltonian

  do iz=0,nz-1	!z
   do icf=0,1	!c/f
    do is=0,1	!up/down 
     do o1=0,l_orb(icf)-1	!orb1

     ind=iz			
     indp=index_zcso(iz,icf,is,o1,nz)

if (icf==1) nf_kr (ind)= nf_kr (ind) + (abs(w_kr_tot(indp,i,ik)))**2 * occup 
if (icf==0) nc_kr (ind)= nc_kr (ind) + (abs(w_kr_tot(indp,i,ik)))**2 * occup	
   
      do iz2=0,nz-1	!z
       do icf2=0,1	!c/f
        do is2=0,1	!up/down 
         do o2=0,l_orb(icf2)-1	!orb2
         
     ind2=iz2			
     ind2p=index_zcso(iz2,icf2,is2,o2,nz)


!same plane
if (icf==1 .and. icf2==1 .and. iz2==iz) nfk_kr (ind)= nfk_kr (ind) +real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
		   	 	                     +conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik)) &
		  	  		   	     *occup/2


if (icf==1 .and. icf2==0 .and. iz2==iz) nfck_kr(ind)=nfck_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
		   	 	                     +conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik) &
		  	  		   	     )*occup

!if (icf==0 .and. icf2==1 .and. iz2==iz) nfck_kr(ind)=nfck_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
!		   	 	                     +conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik) &
!		  	  		   	     )*occup

!different plane
if (icf==1 .and. icf2==1 .and. iz2 .ne. iz) nff_kr(ind)=nff_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
			   	 	                +conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik)) &
		  	  		   		*occup*b_kr(ind2)

if (icf==1 .and. icf2==0 .and. iz2 .ne. iz) nfc_kr(ind)=nfc_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
		   	 	                	+conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik) &
		  	  		   		)*occup

!if (icf==0 .and. icf2==1 .and. iz2 .ne. iz) nfc_kr(ind)=nfc_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr_tot(indp,ind2p,ik)*w_kr_tot(ind2p,i,ik)&
!		   	 	                	+conjg(w_kr_tot(ind2p,i,ik))*V_kr_tot(ind2p,indp,ik)*w_kr_tot(indp,i,ik) &
!		  	  		   		)*occup






	

!if (icf==1 .and. icf2==0 .and. iz2==iz) nfck_kr(ind)=nfck_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr(indp,ind2p)*w_kr_tot(ind2p,i,ik)&
!		   	 	                     +conjg(w_kr_tot(ind2p,i,ik))*V_kr(ind2p,indp)*w_kr_tot(indp,i,ik)) &
!		  	  		   	     *occup*k_vecs2d(ik,4)

!if (icf==1 .and. icf2==0 .and. iz2 .ne. iz) nfc_kr(ind)=nfc_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr(indp,ind2p)*w_kr_tot(ind2p,i,ik)&
!		   	 	                	+conjg(w_kr_tot(ind2p,i,ik))*V_kr(ind2p,indp)*w_kr_tot(indp,i,ik)) &
!		  	  		   		*occup*k_vecs2d(ik,4)

!if (icf==1 .and. icf2==1 .and. iz2 .ne. iz) nff_kr(ind)=nff_kr(ind)+real( conjg(w_kr_tot(indp,i,ik)) *V_kr(indp,ind2p)*w_kr_tot(ind2p,i,ik)&
!			   	 	                +conjg(w_kr_tot(ind2p,i,ik))*V_kr(ind2p,indp)*w_kr_tot(indp,i,ik)) &
!		  	  		   		*occup*b_kr(ind2)*k_vecs2d(ik,4)


!if (icf==1 .and. icf2==1) nff(ind)=nff(ind)+real( conjg(w(indp,i)) *h_psi_tot(indp,j)*w(ind2p,i)&
!		   	 	                 +conjg(w(ind2p,i))*conjg(h_psi_tot(indp,j))*w(indp,i)) &
!		   	 	            *b(ind2) &
!		  	  		   *occup


!!!nfc
        
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo	!eigenvalues
enddo	!k

!fctot=sum(nfc(0:n_sites-1))	!total <Vf^+c + V*c^+f>

Nftot=sum(nf_kr(0:n_sites_kr-1))	!total number of f electrons
Nctot=sum(nc_kr(0:n_sites_kr-1))	!total number of c electrons
Nfctot=Nftot+Nctot		!total number of f+c electrons




nf_kr=nf_kr/n_k_kr
nc_kr=nc_kr/n_k_kr
nff_kr=nff_kr/n_k_kr
nfc_kr=nfc_kr/n_k_kr
nfk_kr=nfk_kr/n_k_kr
nfck_kr=nfck_kr/n_k_kr

!nfctot=(sum(nf_kr(0:n_sites_kr-1))+sum(nc_kr(0:n_sites_kr-1)))*nkx*nky

!!!!!!!!!total free energy
!E_tot_k=0

!from eigenvalues
!do i=0,3
! do j=0,n_k-1
!  E_tot_k=E_tot_k+e_k_tot(i,j)/(exp((e_k_tot(i,j)-mu_new_k)/T)+1)
! enddo
!enddo

!from b
!E_tot_k=E_tot_k-n_k*lambda_new_k*(b_new_k**2-1)

!from mu
!E_tot_k=E_tot_k+mu_new_k*N_el_k

!from lambda**2
!E_tot_k_l=E_tot_k-n_k*lambda_new_k**2/(2*lambda_par_k)

!E_tot_k=E_tot_k/n_k
!E_tot_k_l=E_tot_k_l/n_k

!!!!!!


!!!!!!set new variables
diff=0
do i=0, nz-1

 if(lvac_kr(i)) then
 
  b_new_kr(i)=0
  lambda_new_kr(i)=100
!  ef_kr(i)=ef_site_kr(i,0)-lambda_new_kr(i)

 else


if (scf_type=='sb') then
!b_new(i)=(nfc(i)+nff(i))/(2*lambda(i))				!r
!b_new_k=(nfc_k/(lambda_k-nff_k)/2)				!k

  b_new_kr(i)=(nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(2*lambda_kr(i)-2*nfk_kr(i))
  lambda_new_kr(i)=lambda_kr(i)+lambda_par_kr*(1-b_new_kr(i)**2-nf_kr(i))
!  ef_kr(i)=ef_site_kr(i,0)-lambda_new_kr(i)!-mu
 
  diff=diff+(1-b_new_kr(i)**2-nf_kr(i))**2/nz+(b_new_kr(i)-(nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(2*lambda_kr(i)-2*nfk_kr(i)))**2

else if (scf_type=='ga') then

!b_new(i)=(nfc(i)+nff(i))/(4*lambda(i))*(2.-b(i)**2)**2		r
!b_new_k=nfc_k/2/(-nff_k+2*lambda_k/(2-b_k**2)**2)		k

  b_new_kr(i)=(nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(4*lambda_kr(i)/(2-b_kr(i)**2)**2-2*nfk_kr(i))
  lambda_new_kr(i)=lambda_kr(i)+lambda_par_kr*((1-b_new_kr(i)**2)/(1-b_new_kr(i)**2/2)-nf_kr(i))
!  ef_kr(i)=ef_site_kr(i,0)-lambda_new_kr(i)!-mu
 
  diff=diff+((1-b_new_kr(i)**2)/(1-b_new_kr(i)**2/2)-nf_kr(i))**2/nz+(b_new_kr(i)-(nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(4*lambda_kr(i)/(2-b_new_kr(i)**2)**2-2*nfk_kr(i)))**2


endif


 endif
 
enddo

if (.not. lfixmu_kr) diff=diff+((Nfctot-N_el_kr)/n_k_kr)**2

diff=sqrt(diff)

!!!!!!!!!!Writes new parameters
ef_kr=0
write(*,"(i6,100f20.12)"), i_iter, diff
write (*, "(a5,21(a10))"), "i", "b_old","nf","b_new","bnfc","ef","b^2+nf","nfc","lambda","l_new","mu_new", "nc", "nf+nc", "nff", &
					      "nfk","nfck","E_tot_k"
do i =0, nz-1
write(*,"(i5,100f10.5)")  i, b_kr(i), nf_kr(i), b_new_kr(i), (nfc_kr(i)+nfck_kr(i)+nff_kr(i))/(2*lambda_kr(i)-2*nfk_kr(i)),&
			  ef_kr(i), b_new_kr(i)**2+nf_kr(i), nfc_kr(i), lambda_kr(i), lambda_new_kr(i), mu_new_kr, nc_kr(i), &
				 nf_kr(i)+nc_kr(i), nff_kr(i),nfk_kr(i),nfck_kr(i)!, E_tot_kr
enddo



write(*,"(10(a18))") "filling", "ntot", "nf/n_sites","nc/n_sites","nf","nc"  
write(*,"(100f18.10)") , Nfctot/n_k_kr/n_sites_kr, Nfctot,Nftot/n_k_kr/n_sites_kr, Nctot/n_k_kr/n_sites_kr, nftot, nctot


end subroutine compute_nf_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_ham_kr(Ham_kr)
complex(KIND=idpc) 	:: Ham_kr(0:size_h_kr-1,0:size_h_kr-1)

do, i=0,size_h_kr-1
    write(*,"(100(a1,f5.2,a1,f5.2,a1))") ( "(",real(Ham_kr(i,j)),",",aimag(Ham_kr(i,j)),")", j=0,size_h_kr-1 )
enddo
   write(*, " ")


end subroutine write_ham_kr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_bands_kr(e_kr_tot)
real (kind=idp)		:: e_kr_tot(0:size_h_kr-1,0:nk2d-1)
!complex(KIND=idpc) 	:: w_kr_tot
integer		::i,j,ik



  write(21,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr-1
     do ik= 0, nk2d-1
      write(21,"(3f10.5, 5f15.10)") k_vecs2d(ik,1),k_vecs2d(ik,2), k_vecs2d(ik,3),e_kr_tot(i,ik)
      enddo
      write(21,"(a1)") ""
     enddo




end subroutine write_bands_kr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_wfc_kr(e_kr_tot,w_kr_tot)
complex(KIND=idpc) 	:: w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1)
real (kind=idp)		:: e_kr_tot(0:size_h_kr-1,0:nk2d-1)
integer 			:: i,j,ix,iy,iz,icf,is, indp1, indp2, indp3, indp4
real(KIND=idpc)			:: e_thr

e_thr=0.1


do i=0, size_h_kr-1		
 do ik=0,nk2d-1

   write (23, "(a1,2i6,f20.10)") "#",i,ik,e_kr_tot(i,ik)
! if(abs(e(i))<e_thr) then
!  write (23, "(a3,i6,2f20.10)") "#",i,e(i)
 
   do iz=0,nz-1	!z
    
     indp1=iz + 0*nz + 0*2*nz	!c up
     indp2=iz + 0*nz + 1*2*nz	!c down
     indp3=iz + 1*nz + 0*2*nz	!f up
     indp4=iz + 1*nz + 1*2*nz	!f down
		     
     write (23, "(3i5,9f15.6)") ikx,iky,iz,abs(w_kr_tot(indp1,i,ik))**2+abs(w_kr_tot(indp2,i,ik))**2+abs(w_kr_tot(indp3,i,ik))**2+abs(w_kr_tot(indp4,i,ik))**2,&
     				 real(w_kr_tot(indp1,i,ik)), aimag(w_kr_tot(indp1,i,ik)),real(w_kr_tot(indp2,i,ik)), aimag(w_kr_tot(indp2,i,ik)),&
     				 real(w_kr_tot(indp3,i,ik)), aimag(w_kr_tot(indp3,i,ik)),real(w_kr_tot(indp4,i,ik)), aimag(w_kr_tot(indp4,i,ik))
     					

   enddo
 
  write (23,"(a1)") ""
  write (23,"(a1)") ""

! endif
 enddo
enddo



end subroutine write_wfc_kr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kvecs_pw (k_vecs3d_irr, nk3d_irr,nkx_k,nky_k,nkz,lysameasx,lzsameasx,nky)
real(kind=idp),allocatable, intent(out)	::k_vecs3d_irr(:,:)
integer, intent(out)	::nk3d_irr
integer	::ik,ikn,kshifted
real(kind=idp)	:: y_corr, z_corr
integer,intent(in):: nkx_k,nky_k,nkz,nky
logical,intent(in):: lysameasx,lzsameasx


if (lysameasx) then
y_corr=1.0d0
else
y_corr=0.9999d0
endif
if (lzsameasx) then
z_corr=1.0d0
else
z_corr=0.9998d0
endif


if (lkshifted) kshifted=1
if (.not. lkshifted) kshifted=0


OPEN(unit=100,file=trim(label)//'/input_k3d',status='unknown')

write(100,"(i1)") 8
write(100,"(a)") "k3d"
write(100,"(f15.8)") real(nkz)/real(nkx_k)*z_corr	!to impose that z direction is always different from x-y
write(100,"(f15.8)") real(nky)/real(nkx_k)*y_corr
write(100,"(3i6)") nkx_k,nky_k,nkz
write(100,"(3i6)") kshifted,kshifted,0
write(100,"(a1)") "f"

 CLOSE(unit=100)


 call system("cat " // trim(label) // "/input_k3d")
 call system("kpoints.x <" // trim(label) // "/input_k3d")
 call system("mv info " // trim(label) )
 call system("mv k3d " // trim(label) )


OPEN(unit=101,file=trim(label)//'/k3d',status='unknown')
read(101,"(i)"), nk3d_irr

allocate(k_vecs3d_irr(0:nk3d_irr-1,4))

do ik=0, nk3d_irr-1
!read(101,"(i5,1x,3f20.10,f7.2)"), ikn, k_vecs3d_irr(ikn-1,1),k_vecs3d_irr(ikn-1,2),k_vecs3d_irr(ikn-1,3),k_vecs3d_irr(ikn-1,4)
read(101,*), ikn, k_vecs3d_irr(ikn-1,1),k_vecs3d_irr(ikn-1,2),k_vecs3d_irr(ikn-1,3),k_vecs3d_irr(ikn-1,4)
enddo
 CLOSE(unit=101)

write(*,"(a21,i5)") "Nk points irrid:", nk3d_irr

do ik=0, nk3d_irr-1
k_vecs3d_irr(ik,2)=k_vecs3d_irr(ik,2)*real(nky)/real(nkx_k)*y_corr
k_vecs3d_irr(ik,3)=k_vecs3d_irr(ik,3)*real(nkz)/real(nkx_k)*z_corr
!write(*, "(i5,1x,3f20.10,f7.2)"), ik, k_vecs3d_irr(ik,1),k_vecs3d_irr(ik,2),k_vecs3d_irr(ik,3),k_vecs3d_irr(ik,4)
enddo


end subroutine set_kvecs_pw


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_k_vecs_nscf!(k_vecs_nscf,factor_k_nscf,phi_k_sigma_alpha_nscf,phi_k_alpha_sigma_nscf, delta_k_nscf)
!real(kind=idp),allocatable,intent(out)	::k_vecs_nscf(:,:),factor_k_nscf(:),delta_k_nscf(:)
!complex(kind=idpc),allocatable,intent(out)	::phi_k_sigma_alpha_nscf(:,:,:,:,:),phi_k_alpha_sigma_nscf(:,:,:,:,:)

if (.not. lpathnscf)  nk_plot=nkx_plot*nky_plot*nkz_plot
if (.not. lpathnscf)  nk_plot=nkx_plot**3
if (.not. lpathnscf .and. lcheckcubic)  nk_plot=13
if (lpathnscf) nk_plane=4*nkpath_plot+2*nkpath_plot_xy+1
if (lpathnscf) nk_plot=nk_plane*nkz_plot


if (.not. lpathnscf .and. .not. lcheckcubic) then
print *, "nopath nocheck"


!lk3d=1
print *,lk3d

if (lk3d==0) then
!call set_kvecs_pw(k_vecs_nscf,nk_plot,nk_dos,nk_dos,nk_dos,.true.,.true.,nkx_plot)	!cubic
 call set_kvecs_pw(k_vecs_nscf,nk_plot,nk_dos,nk_dos,nk_dos,.true.,.false.,nkx_plot)	!tetragonal

total_weight=sum(k_vecs_nscf(:,4))
k_vecs_nscf(:,4)=k_vecs_nscf(:,4)/total_weight

elseif (lk3d==1) then
OPEN(unit=101,file='k3d_50',status='unknown')
read(101,"(i)"), nk_plot

allocate(k_vecs_nscf(0:nk_plot-1,4))

do ik=0, nk_plot-1
read(101,"(i5,1x,3f20.10,f7.2)"), ikn, k_vecs_nscf(ikn-1,1),k_vecs_nscf(ikn-1,2),k_vecs_nscf(ikn-1,3),k_vecs_nscf(ikn-1,4)
enddo
 CLOSE(unit=101)

!adds normalization
total_weight=0
do ik=0, nk_plot-1
total_weight=total_weight+k_vecs_nscf(ik,4)
enddo

do ik=0, nk_plot-1
k_vecs_nscf(ik,4)=k_vecs_nscf(ik,4)/total_weight
enddo



elseif (lk3d==2) then

 call set_kvecs_cubic(k_vecs_nscf,nk_plot,nk_dos)	!cubic

elseif (lk3d==3) then		!ok for all symmetries (points in the whole BZ)

nk_plot=nk_dos**3
allocate(k_vecs_nscf(0:nk_plot-1,4))

do ikx=0,nk_dos-1
do iky=0,nk_dos-1
  do ikz=0,nk_dos-1
   	ik=ikx + nk_dos*iky + ikz*nk_dos**2
	
if (mod(nk_dos,2)==1)   k_vecs_nscf (ik,1)=-0.5+real(ikx+0.5)/real(nk_dos)	!k_x 	odd
if (mod(nk_dos,2)==1)   k_vecs_nscf (ik,2)=-0.5+real(iky+0.5)/real(nk_dos)	!k_y
if (mod(nk_dos,2)==1)   k_vecs_nscf (ik,3)=-0.5+real(ikz+0.5)/real(nk_dos)	!k_z

if (mod(nk_dos,2)==0)   k_vecs_nscf (ik,1)=-0.5+real(ikx)/real(nk_dos)	!k_x	even
if (mod(nk_dos,2)==0)   k_vecs_nscf (ik,2)=-0.5+real(iky)/real(nk_dos)	!k_y
if (mod(nk_dos,2)==0)   k_vecs_nscf (ik,3)=-0.5+real(ikz)/real(nk_dos)	!k_z

k_vecs_nscf(ik,4)=1.0d0/dble(nk_plot)

  enddo
 enddo
enddo


elseif (lk3d==4) then		!ok for all symmetries (points in the whole BZ) tetragonal sampling to mimic the kr procedure

nk_plot=nk_dos**2*nz
allocate(k_vecs_nscf(0:nk_plot-1,4))

do ikx=0,nk_dos-1
do iky=0,nk_dos-1
  do ikz=0,nz-1
   	ik=ikx + nk_dos*iky + ikz*nk_dos**2
	
if (mod(nk_dos,2)==1)  	 k_vecs_nscf (ik,1)=-0.5+real(ikx+0.5)/real(nk_dos)	!k_x 	odd
if (mod(nk_dos,2)==1)  	 k_vecs_nscf (ik,2)=-0.5+real(iky+0.5)/real(nk_dos)	!k_y
if (mod(nz,2)==1)  	 k_vecs_nscf (ik,3)=-0.5+real(ikz+0.5)/real(nz)	!k_z

if (mod(nk_dos,2)==0)  	 k_vecs_nscf (ik,1)=-0.5+real(ikx)/real(nk_dos)	!k_x	even
if (mod(nk_dos,2)==0)  	 k_vecs_nscf (ik,2)=-0.5+real(iky)/real(nk_dos)	!k_y
if (mod(nz,2)==0) 	 k_vecs_nscf (ik,3)=-0.5+real(ikz)/real(nz)	!k_z

k_vecs_nscf(ik,4)=1.0d0/dble(nk_plot)

  enddo
 enddo
enddo



endif


elseif (.not. lpathnscf .and. lcheckcubic) then
allocate(k_vecs_nscf(0:nk_plot-1,4))
print *, "nopath check"
k1_c=0.1
k2_c=0.2
k3_c=0.3

k_vecs_nscf=0
!permutations
k_vecs_nscf(0,1)=k1_c
k_vecs_nscf(0,2)=k2_c
k_vecs_nscf(0,3)=k3_c
k_vecs_nscf(1,1)=k1_c
k_vecs_nscf(1,2)=k3_c
k_vecs_nscf(1,3)=k2_c
k_vecs_nscf(2,1)=k2_c
k_vecs_nscf(2,2)=k1_c
k_vecs_nscf(2,3)=k3_c
k_vecs_nscf(3,1)=k2_c
k_vecs_nscf(3,2)=k3_c
k_vecs_nscf(3,3)=k1_c
k_vecs_nscf(4,1)=k3_c
k_vecs_nscf(4,2)=k1_c
k_vecs_nscf(4,3)=k2_c
k_vecs_nscf(5,1)=k3_c
k_vecs_nscf(5,2)=k2_c
k_vecs_nscf(5,3)=k1_c
!reflections
k_vecs_nscf(6,1)=-k1_c
k_vecs_nscf(6,2)=k2_c
k_vecs_nscf(6,3)=k3_c
k_vecs_nscf(7,1)=k1_c
k_vecs_nscf(7,2)=-k2_c
k_vecs_nscf(7,3)=k3_c
k_vecs_nscf(8,1)=k1_c
k_vecs_nscf(8,2)=k2_c
k_vecs_nscf(8,3)=-k3_c
k_vecs_nscf(9,1)=-k1_c
k_vecs_nscf(9,2)=-k2_c
k_vecs_nscf(9,3)=k3_c
k_vecs_nscf(10,1)=-k1_c
k_vecs_nscf(10,2)=k2_c
k_vecs_nscf(10,3)=-k3_c
k_vecs_nscf(11,1)=k1_c
k_vecs_nscf(11,2)=-k2_c
k_vecs_nscf(11,3)=-k3_c
k_vecs_nscf(12,1)=-k1_c
k_vecs_nscf(12,2)=-k2_c
k_vecs_nscf(12,3)=-k3_c
!endif


allocate(phi_k_sigma_alpha_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1))
allocate(phi_k_alpha_sigma_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot-1))
allocate(kin_energy_k_nscf(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1,0:1))
allocate(delta_k_nscf(0:nk_plot-1))
allocate(fermien_k(0:nk_plot-1,2))
allocate(occup_k(0:nk_plot-1))



else		!path
print *, "path"

allocate(k_vecs_nscf(0:nk_plot-1,4))
allocate(phi_k_sigma_alpha_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1))
allocate(phi_k_alpha_sigma_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot-1))
allocate(kin_energy_k_nscf(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1,0:1))
allocate(delta_k_nscf(0:nk_plot-1))
allocate(fermien_k(0:nk_plot-1,2))
allocate(occup_k(0:nk_plot-1))


do ikz=0, nkz_plot-1

 do ikx=0, nk_plane-1
  ik= ikx + nk_plane*ikz
  if (mod(nkz_plot,2)==1) k_vecs_nscf (ik,3)=-0.5+real(ikz+0.5)/real(nkz_plot)	!k_z
  if (mod(nkz_plot,2)==0) k_vecs_nscf (ik,3)=-0.5+real(ikz)/real(nkz_plot)	!k_z
 enddo


  do ikx=0,nkpath_plot-1

   ik= ikx + nk_plane*ikz
  
   k_vecs_nscf (ik,1)=real(ikx)/real(2*nkpath_plot)
   k_vecs_nscf (ik,2)=0

   k_vecs_nscf (nkpath_plot+ik,1)=0.5
   k_vecs_nscf (nkpath_plot+ik,2)=real(ikx)/real(2*nkpath_plot)

   k_vecs_nscf (2*nkpath_plot+nkpath_plot_xy+ik,1)=0
   k_vecs_nscf (2*nkpath_plot+nkpath_plot_xy+ik,2)=real(ikx)/real(2*nkpath_plot)

   k_vecs_nscf (3*nkpath_plot+2*nkpath_plot_xy+ik,1)=-0.5+real(ikx)/real(2*nkpath_plot)
   k_vecs_nscf (3*nkpath_plot+2*nkpath_plot_xy+ik,2)=0

  enddo

 do ikx=0,nkpath_plot_xy-1
  
   ik= ikx + nk_plane*ikz
   
   k_vecs_nscf (2*nkpath_plot+ik,1)=0.5-real(ikx)/real(2*nkpath_plot_xy)
   k_vecs_nscf (2*nkpath_plot+ik,2)=k_vecs_nscf (2*nkpath_plot+ik,1)

   k_vecs_nscf (3*nkpath_plot+nkpath_plot_xy+ik,1)=-real(ikx)/real(2*nkpath_plot_xy)
   k_vecs_nscf (3*nkpath_plot+nkpath_plot_xy+ik,2)=0.5-real(ikx)/real(2*nkpath_plot_xy)

 enddo

! k_vecs_nscf (3*nkpath_plot+2*nkpath_plot_xy+1,1)=0
! k_vecs_nscf (3*nkpath_plot+2*nkpath_plot_xy+1,2)=0

enddo

endif

if (.not. allocated(phi_k_sigma_alpha_nscf)) allocate(phi_k_sigma_alpha_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1))
if (.not. allocated(phi_k_alpha_sigma_nscf)) allocate(phi_k_alpha_sigma_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot-1))
if (.not. allocated(kin_energy_k_nscf)) allocate(kin_energy_k_nscf(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot-1,0:1))
if (.not. allocated(delta_k_nscf)) allocate(delta_k_nscf(0:nk_plot-1))
if (.not. allocated(fermien_k)) allocate(fermien_k(0:nk_plot-1,2))
if (.not. allocated(occup_k)) allocate(occup_k(0:nk_plot-1))



!call factor_K(factor_k_nscf, nk_plot, k_vecs_nscf,3)

call set_phi_k(phi_k_sigma_alpha_nscf,phi_k_alpha_sigma_nscf, kin_energy_k_nscf,delta_k_nscf,nk_plot,k_vecs_nscf,3)


!do ik=0,nk_plot-1
!  write(*,"(3f10.5, 100f7.4)") k_vecs_nscf(ik,1),k_vecs_nscf(ik,2),k_vecs_nscf(ik,3), factor_k_nscf(ik),&
!  				real(phi_k_sigma_alpha_nscf(0,0,ik)),real(phi_k_sigma_alpha_nscf(0,1,ik)),real(phi_k_sigma_alpha_nscf(1,0,ik)),real(phi_k_sigma_alpha_nscf(1,1,ik)),&
!  				aimag(phi_k_sigma_alpha_nscf(0,0,ik)),aimag(phi_k_sigma_alpha_nscf(0,1,ik)),aimag(phi_k_sigma_alpha_nscf(1,0,ik)),aimag(phi_k_sigma_alpha_nscf(1,1,ik))
!enddo



end subroutine set_k_vecs_nscf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_k_vecs_chern(iplane)
integer::iplane


if (iplane==1 .or. iplane==2) nk_chern=nkx_chern*nky_chern	!k_z=0 or k_z=pi
if (iplane==3) then						!k_x=k_y
nkxy_chern=int(nkx_chern*sqrt(2.0d0))
nkz_chern=nkx_chern
nk_chern=nkxy_chern*nkz_chern
endif

allocate(k_vecs_chern(0:nk_chern-1,4))
allocate(phi_k_sigma_alpha_chern(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_chern-1))
allocate(phi_k_alpha_sigma_chern(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_chern-1))
allocate(kin_energy_k_chern(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nk_chern-1,0:1))
allocate(delta_k_chern(0:nk_chern-1))


k_vecs_chern=0


if(iplane==1 .or. iplane==2) then

if (iplane_xyz==1) then
i=2	!y
j=3	!z
k=1	!x=0,pi
elseif(iplane_xyz==2) then
i=3	!z	
j=1	!x
k=2	!y=0,pi
elseif (iplane_xyz==3) then
i=1	!x
j=2	!y
k=3	!z=0,pi
endif

do ikx=0,nkx_chern-1
 do iky=0,nky_chern-1
   	ik=ikx + nkx_chern*iky
	
if (mod(nkx_chern,2)==0)   k_vecs_chern (ik,i)=-0.5+real(ikx)/real(nkx_chern)	!k_x 	odd
if (mod(nky_chern,2)==0)   k_vecs_chern (ik,j)=-0.5+real(iky)/real(nky_chern)	!k_y

if (mod(nkx_chern,2)==1)   k_vecs_chern (ik,i)=-0.5+real(ikx+0.5)/real(nkx_chern)	!k_x	even
if (mod(nky_chern,2)==1)   k_vecs_chern (ik,j)=-0.5+real(iky+0.5)/real(nky_chern)	!k_y

if (iplane==1) k_vecs_chern (ik,k)=0.0d0	!k_z
if (iplane==2) k_vecs_chern (ik,k)=0.5d0	!k_z

  enddo
 enddo


elseif (iplane==3) then

!if (iplane_xyz==3) then

do ikx=0,nkxy_chern-1
 do ikz=0,nkz_chern-1
   	ik=ikx + nkxy_chern*ikz
	
if (mod(nkxy_chern,2)==0)   k_vecs_chern (ik,1)=-0.5+real(ikx)/real(nkxy_chern)	!k_x 	odd
!if (mod(nkxy_chern,2)==0)   k_vecs_chern (ik,1)=-0.5+real(ikx+0.5)/real(nkxy_chern)	!k_x 	odd
if (mod(nkz_chern,2)==0)    k_vecs_chern (ik,3)=-0.5+real(ikz)/real(nkz_chern)	!k_y

if (mod(nkxy_chern,2)==1)   k_vecs_chern (ik,1)=-0.5+real(ikx+0.5)/real(nkxy_chern)	!k_x	even
!if (mod(nkxy_chern,2)==1)   k_vecs_chern (ik,1)=-0.5+real(ikx)/real(nkxy_chern)	!k_x	even
if (mod(nkz_chern,2)==1)    k_vecs_chern (ik,3)=-0.5+real(ikz+0.5)/real(nkz_chern)	!k_y

k_vecs_chern (ik,2)=k_vecs_chern (ik,1)

  enddo
 enddo

!elseif (iplane_xyz==1) then
!do ikx=0,nkxy_chern-1
! do ikz=0,nkz_chern-1
!   	ik=ikx + nkxy_chern*ikz
	
!if (mod(nkxy_chern,2)==0)   k_vecs_chern (ik,2)=-0.5+real(ikx)/real(nkxy_chern)	!k_x 	odd
!if (mod(nkz_chern,2)==0)    k_vecs_chern (ik,1)=-0.5+real(ikz)/real(nkz_chern)	!k_y

!if (mod(nkxy_chern,2)==1)   k_vecs_chern (ik,2)=-0.5+real(ikx+0.5)/real(nkxy_chern)	!k_x	even
!if (mod(nkz_chern,2)==1)    k_vecs_chern (ik,1)=-0.5+real(ikz+0.5)/real(nkz_chern)	!k_y

!k_vecs_chern (ik,3)=k_vecs_chern (ik,2)

!  enddo
! enddo

!elseif (iplane_xyz==2) then
!do ikx=0,nkxy_chern-1
! do ikz=0,nkz_chern-1
!   	ik=ikx + nkxy_chern*ikz
	
!if (mod(nkxy_chern,2)==0)   k_vecs_chern (ik,3)=-0.5+real(ikx)/real(nkxy_chern)	!k_x 	odd
!if (mod(nkz_chern,2)==0)    k_vecs_chern (ik,2)=-0.5+real(ikz)/real(nkz_chern)	!k_y

!if (mod(nkxy_chern,2)==1)   k_vecs_chern (ik,3)=-0.5+real(ikx+0.5)/real(nkxy_chern)	!k_x	even
!if (mod(nkz_chern,2)==1)    k_vecs_chern (ik,2)=-0.5+real(ikz+0.5)/real(nkz_chern)	!k_y

!k_vecs_chern (ik,1)=k_vecs_chern (ik,3)

!  enddo
! enddo
!endif

!do ik=0, nk_chern-1
!write (*, "(i8, 20f10.5)"), ik,  k_vecs_chern (ik,1), k_vecs_chern (ik,2), k_vecs_chern (ik,3)
!enddo
endif

!!!!!!!!!!!!!!then sets ham parts

call set_phi_k(phi_k_sigma_alpha_chern,phi_k_alpha_sigma_chern, kin_energy_k_chern,delta_k_chern,nk_chern,k_vecs_chern,3)


end subroutine set_k_vecs_chern

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_k_vecs_path3d!(k_vecs_path3d,factor_k_path3d,phi_k_sigma_alpha_path3d,phi_k_alpha_sigma_path3d, delta_k_path3d)
!real(kind=idp),allocatable,intent(out)	::k_vecs_path3d(:,:),factor_k_path3d(:),delta_k_path3d(:)
!complex(kind=idpc),allocatable,intent(out)	::phi_k_sigma_alpha_path3d(:,:,:),phi_k_alpha_sigma_path3d(:,:,:)
integer	::nkpath_plot_xyz

nkpath_plot_xyz=int(sqrt(3.0D0)*nkpath_plot)

if(1) then

nk_plot_path3d=4*nkpath_plot+2*nkpath_plot_xy+nkpath_plot_xyz+1

print *, "nk_plot_path3d", nk_plot_path3d

allocate(k_vecs_path3d(0:nk_plot_path3d-1,4))
!allocate(factor_k_path3d(0:nk_plot_path3d-1))
allocate(phi_k_sigma_alpha_path3d(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot_path3d-1))
allocate(phi_k_alpha_sigma_path3d(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot_path3d-1))
allocate(kin_energy_k_path3d(0:l_orbital_f-1,0:l_orbital_f-1, 0:1,0:1,0:nk_plot_path3d-1,0:1))
allocate(delta_k_path3d(0:nk_plot_path3d-1))
!allocate(fermien_k(0:nk_plot-1,2))
!allocate(occup_k(0:nk_plot-1))


k_vecs_path3d=0

do ikx=0,nkpath_plot-1
 
!(0,0,0)->(0,0,pi) GX
   k_vecs_path3d (ikx,1)=0
   k_vecs_path3d (ikx,2)=0
   k_vecs_path3d (ikx,3)=dble(ikx)/dble(2*nkpath_plot)
!   k_vecs_path3d (ikx,1)=0
!   k_vecs_path3d (ikx,2)=dble(ikx)/dble(2*nkpath_plot)
!   k_vecs_path3d (ikx,3)=0

!(0,0,pi)->(pi,0,pi) XM
   k_vecs_path3d (nkpath_plot+ikx,1)=dble(ikx)/dble(2*nkpath_plot)
   k_vecs_path3d (nkpath_plot+ikx,2)=0
   k_vecs_path3d (nkpath_plot+ikx,3)=0.5d0

enddo

do ikx=0,nkpath_plot_xy-1
  
!(pi,0,pi)->(0,0,0) MG
   k_vecs_path3d (2*nkpath_plot+ikx,1)=0.5d0-dble(ikx)/dble(2*nkpath_plot_xy)
   k_vecs_path3d (2*nkpath_plot+ikx,2)=0
   k_vecs_path3d (2*nkpath_plot+ikx,3)=0.5d0-dble(ikx)/dble(2*nkpath_plot_xy)

enddo

do ikx=0,nkpath_plot_xyz-1
 
!(0,0,0)->(pi,pi,pi) GR
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+ikx,1)=dble(ikx)/dble(2*nkpath_plot_xyz)
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+ikx,2)=dble(ikx)/dble(2*nkpath_plot_xyz)
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+ikx,3)=dble(ikx)/dble(2*nkpath_plot_xyz)

enddo

do ikx=0,nkpath_plot-1
 
!(pi,pi,pi)->(0,pi,pi) RM
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,1)=0.5d0-dble(ikx)/dble(2*nkpath_plot)
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,2)=0.5d0
   k_vecs_path3d (2*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,3)=0.5d0

!(0,pi,pi)->(0,0,pi) MX
   k_vecs_path3d (3*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,1)=0
   k_vecs_path3d (3*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,2)=0.5d0-dble(ikx)/dble(2*nkpath_plot)
   k_vecs_path3d (3*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,3)=0.5d0

enddo


do ikx=0,nkpath_plot_xy-1
  
!(0,0,pi)->(pi,pi,pi) XR
   k_vecs_path3d (4*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,1)=dble(ikx)/dble(2*nkpath_plot_xy)
   k_vecs_path3d (4*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,2)=dble(ikx)/dble(2*nkpath_plot_xy)
   k_vecs_path3d (4*nkpath_plot+nkpath_plot_xy+nkpath_plot_xyz+ikx,3)=0.5d0

enddo

   k_vecs_path3d (4*nkpath_plot+2*nkpath_plot_xy+nkpath_plot_xyz,1)=0.5d0
   k_vecs_path3d (4*nkpath_plot+2*nkpath_plot_xy+nkpath_plot_xyz,2)=0.5d0
   k_vecs_path3d (4*nkpath_plot+2*nkpath_plot_xy+nkpath_plot_xyz,3)=0.5d0


else

nk_plot_path3d=3*nkpath_plot+nkpath_plot_xyz+1

print *, "nk_plot_path3d", nk_plot_path3d

allocate(k_vecs_path3d(0:nk_plot_path3d-1,3))
!allocate(factor_k_path3d(0:nk_plot_path3d-1))
allocate(phi_k_sigma_alpha_path3d(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot_path3d-1))
allocate(phi_k_alpha_sigma_path3d(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot_path3d-1))
allocate(kin_energy_k_path3d(0:l_orbital_f-1,0:l_orbital_f-1, 0:1,0:1,0:nk_plot_path3d-1,0:1))
allocate(delta_k_path3d(0:nk_plot_path3d-1))
!allocate(fermien_k(0:nk_plot-1,2))
!allocate(occup_k(0:nk_plot-1))


k_vecs_path3d=0

do ikx=0,nkpath_plot-1
 
!(0,0,0)->(0,0,pi)
   k_vecs_path3d (ikx,1)=0
   k_vecs_path3d (ikx,2)=0
   k_vecs_path3d (ikx,3)=dble(ikx)/dble(2*nkpath_plot)

!(0,0,pi)->(0,pi,pi)
   k_vecs_path3d (nkpath_plot+ikx,1)=0
   k_vecs_path3d (nkpath_plot+ikx,2)=dble(ikx)/dble(2*nkpath_plot)
   k_vecs_path3d (nkpath_plot+ikx,3)=0.5d0


!(0,pi,pi)->(pi,pi,pi)
   k_vecs_path3d (2*nkpath_plot+ikx,1)=dble(ikx)/dble(2*nkpath_plot)
   k_vecs_path3d (2*nkpath_plot+ikx,2)=0.5d0
   k_vecs_path3d (2*nkpath_plot+ikx,3)=0.5d0

enddo

!(pi,pi,pi)->(0,0,0)

do ikx=0,nkpath_plot_xyz
 
   k_vecs_path3d (3*nkpath_plot+ikx,1)=0.5d0-dble(ikx)/dble(2*nkpath_plot_xyz)
   k_vecs_path3d (3*nkpath_plot+ikx,2)=0.5d0-dble(ikx)/dble(2*nkpath_plot_xyz)
   k_vecs_path3d (3*nkpath_plot+ikx,3)=0.5d0-dble(ikx)/dble(2*nkpath_plot_xyz)

enddo



endif


 
!call factor_K(factor_k_path3d, nk_plot_path3d, 3, k_vecs_path3d,3)

call set_phi_k(phi_k_sigma_alpha_path3d,phi_k_alpha_sigma_path3d, kin_energy_k_path3d, delta_k_path3d,nk_plot_path3d,k_vecs_path3d,3)


!do ik=0,nk_plot_path3d-1
!  write(*,"(i5,3f10.5, 100f7.4)") ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2),k_vecs_path3d(ik,3), factor_k_path3d(ik),&
!  				real(phi_k_sigma_alpha_path3d(0,0,ik)),real(phi_k_sigma_alpha_path3d(0,1,ik)),real(phi_k_sigma_alpha_path3d(1,0,ik)),real(phi_k_sigma_alpha_path3d(1,1,ik)),&
!  				aimag(phi_k_sigma_alpha_path3d(0,0,ik)),aimag(phi_k_sigma_alpha_path3d(0,1,ik)),aimag(phi_k_sigma_alpha_path3d(1,0,ik)),aimag(phi_k_sigma_alpha_path3d(1,1,ik))
!enddo



end subroutine set_k_vecs_path3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_k_vecs_hsp


nk_plot_hsp=8
allocate(k_vecs_hsp(0:nk_plot_hsp-1,4))
allocate(phi_k_sigma_alpha_hsp(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nk_plot_hsp-1))
allocate(phi_k_alpha_sigma_hsp(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nk_plot_hsp-1))
allocate(kin_energy_k_hsp(0:l_orbital_f-1,0:l_orbital_f-1, 0:1,0:1,0:nk_plot_hsp-1,0:1))
allocate(delta_k_hsp(0:nk_plot_hsp-1))
!allocate(fermien_k(0:nk_plot-1,2))
!allocate(occup_k(0:nk_plot-1))


k_vecs_hsp=0
k_vecs_hsp(:,4)=1

!(0,0,0)
k_vecs_hsp(0,1)=0
k_vecs_hsp(0,2)=0
k_vecs_hsp(0,3)=0
!(pi,0,0)
k_vecs_hsp(1,1)=0.5d0
k_vecs_hsp(1,2)=0
k_vecs_hsp(1,3)=0
!(0,pi,0)
k_vecs_hsp(2,1)=0
k_vecs_hsp(2,2)=0.5d0
k_vecs_hsp(2,3)=0
!(0,0,pi)
k_vecs_hsp(3,1)=0
k_vecs_hsp(3,2)=0
k_vecs_hsp(3,3)=0.5d0
!(pi,pi,0)
k_vecs_hsp(4,1)=0.5d0
k_vecs_hsp(4,2)=0.5d0
k_vecs_hsp(4,3)=0
!(pi,0,pi)
k_vecs_hsp(5,1)=0.5d0
k_vecs_hsp(5,2)=0
k_vecs_hsp(5,3)=0.5d0
!(0,pi,pi)
k_vecs_hsp(6,1)=0
k_vecs_hsp(6,2)=0.5d0
k_vecs_hsp(6,3)=0.5d0
!(pi,pi,pi)
k_vecs_hsp(7,1)=0.5d0
k_vecs_hsp(7,2)=0.5d0
k_vecs_hsp(7,3)=0.5d0


call set_phi_k(phi_k_sigma_alpha_hsp,phi_k_alpha_sigma_hsp, kin_energy_k_hsp, delta_k_hsp,nk_plot_hsp,k_vecs_hsp,3)



end subroutine set_k_vecs_hsp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine do_k_nscf(b_k,lambda_k, mu_k)
real(kind=idp),intent(in)	:: b_k,lambda_k, mu_k
integer				:: i,ind
logical				::l_1l


  write(*,"(a30)") "Doing NSCF with:"
  write(*,"(a10,f15.8)") "b=",b_k
  write(*,"(a10,f15.8)") "lambda=",lambda_k
  write(*,"(a10,f15.8)") "mu=",mu_k
  write(*,"(a10,i15)") "k points=",nk_plot


  allocate(e_k_tot_nscf(0:dim_hk-1,0:nk_plot-1))
  allocate(ew_k_tot_nscf(0:dim_hk-1,0:nk_plot-1))
  allocate(w_k_tot_nscf(0:dim_hk-1,0:dim_hk-1,0:nk_plot-1))
  allocate(e_k_tot_path3d(0:dim_hk-1,0:nk_plot_path3d-1))
  allocate(ew_k_tot_path3d(0:dim_hk-1,0:nk_plot_path3d-1,0:l_orb_tot-1))
!  allocate(corr_fct_k_nscf(0:3,0:3))
  allocate(V_k_nscf(0:dim_hk-1,0:dim_hk-1))
!  allocate(e_corr_fct_k_nscf(0:3,0:nk_plot-1))
  
  !write(26,"(3a10, 5a15)") "#k_x","k_y","k_z", "e1", "e2", "e3", "e4"

!fermien_k(:,1)=100		!gives the energy closest to the fermi energy at a given k point
!fermien_k(:,2)=-100		!gives the energy closest to the fermi energy at a given k point
!occup_k=0
ew_k_tot_nscf=0


do ik=0, nk_plot-1
if (mod(ik,10000)==0 .or. ik==nk_plot-1) print *, "k=", ik

!!!builds and diag hamiltonian

  call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk_plot,phi_k_sigma_alpha_nscf,phi_k_alpha_sigma_nscf, kin_energy_k_nscf, k_vecs_nscf)
  
!  call write_ham_k(ham_k)
  
  call diagonalize_ham_k(Ham_k,e_k,w_k)


!  e_k=e_k-mu_k
  e_k_tot_nscf(:,ik)=e_k(:)
  w_k_tot_nscf(:,:,ik)=w_k(:,:)

ew_k_tot_nscf(:,ik)=0

do i=0,2*l_orbital-1
ew_k_tot_nscf(:,ik)=ew_k_tot_nscf(:,ik)+abs(w_k(i,:))**2	!c weight
enddo


enddo

write(24,"(a20)") "K points nscf"
do ik=0, min(nk_plot,100000) -1 

 write(24,"(i6,6f10.5)")  ik, 	 k_vecs_nscf(ik,1), k_vecs_nscf(ik,2), k_vecs_nscf(ik,3), k_vecs_nscf(ik,4)


     do is=0,1
!     write(24,"(4 ('(',f6.2,',',f6.2,')','  ') )")  real(phi_k_sigma_alpha_nscf(is,0,ik)),aimag(phi_k_sigma_alpha_nscf(is,0,ik)), &
 !    						   real(phi_k_sigma_alpha_nscf(is,1,ik)),aimag(phi_k_sigma_alpha_nscf(is,1,ik)),&
  ! 						   real(phi_k_alpha_sigma_nscf(is,0,ik)),aimag(phi_k_alpha_sigma_nscf(is,0,ik)), &
   !						   real(phi_k_alpha_sigma_nscf(is,1,ik)),aimag(phi_k_alpha_sigma_nscf(is,1,ik))
    enddo 

!  write(24, " ")
enddo
write(24,"(a10)") "**********"



  write(26,"(3a10, 10a15)") "#k_x","k_y","k_z", "E", "c weight", "Delta"

!if (.not. lpathnscf) then
!    do i=0,dim_hk-1
!     do iky= 0, nky_plot-1
!      do ikz= 0, nkz_plot-1
!       do ikx= 0, nkx_plot-1
!        ik=ikx+nkx_plot*iky+nkx_plot*nky_plot*ikz
!        write(35,"(3f10.5, 5f15.10)") k_vecs_nscf(ik,1),k_vecs_nscf(ik,2), k_vecs_nscf(ik,3),e_k_tot_nscf(i,ik),ew_k_tot_nscf(i,ik),delta_k_nscf(ik)!,e_corr_fct_k_nscf(i,ik)
!       enddo
!       write(26,"(a10)") "   "
!      enddo
!  write(26,"(a10)") "   "
!     enddo
!    enddo

!     do iky= 0, nky_plot-1
!     do ikz= 0, nkz_plot-1
!     do ikx= 0, nkx_plot-1
!ik=ikx+nkx_plot*iky+nkx_plot*nky_plot*ikz
!if (.not. lpathnscf)      write(33,"(3f10.5, 5f30.15)") k_vecs_nscf(ik,1),k_vecs_nscf(ik,2), k_vecs_nscf(ik,3),&
!			  				min(abs(fermien_k(ik,1)),abs(fermien_k(ik,2))),     fermien_k(ik,1),fermien_k(ik,2), occup_k(ik)
!     enddo
!write(33,"(a1)") " "
!     enddo
!     enddo


!else
if (lpathnscf) then

do i=0,dim_hk-1
 do ikz=0, nkz_plot-1
  do ikx=0,nk_plane-1
 
  ik= ikx + nk_plane*ikz
      write(26,"(i5,3f10.5, 5f15.10)") ikx, k_vecs_nscf(ik,1),k_vecs_nscf(ik,2), k_vecs_nscf(ik,3),e_k_tot_nscf(i,ik),ew_k_tot_nscf(i,ik),delta_k_nscf(ik)!,e_corr_fct_k_nscf(i,ik)

   enddo
  write(26,"(a10)") "   "
  write(26,"(a10)") "   "
  enddo
 enddo
endif

if (lcheckcubic .and. .not. lpathnscf) then

OPEN(unit=500,file=trim(label)//'/check_cubic',status='unknown')

do i=0,dim_hk-1
do ik=0, nk_plot-1 	!k_point
 
      write(500,"(i5,3f10.5, f25.18)") ik, k_vecs_nscf(ik,1),k_vecs_nscf(ik,2), k_vecs_nscf(ik,3),e_k_tot_nscf(i,ik)!,ew_k_tot_nscf(i,ik),delta_k_nscf(ik)!,e_corr_fct_k_nscf(i,ik)

   enddo
  write(500,"(a10)") "   "
  write(500,"(a10)") "   "
 enddo
endif

 close(unit=500)

!DOS
if (.not. lpathnscf .and. .not. lcheckcubic) then
print *, "DOS"
e_min=omega_min
e_max=omega_max
delta_e=delta_omega

OPEN(unit=500,file=trim(label)//'/dos_k_3d',status='unknown')
OPEN(unit=600,file=trim(label)//'/g_k_3d',status='unknown')
OPEN(unit=800,file=trim(label)//'/g_k_3d_all',status='unknown')
OPEN(unit=700,file=trim(label)//'/g_k_3d_kres',status='unknown')
allocate(dos_k(0:dim_hk+3))
allocate(g_k(0:dim_hk-1,0:dim_hk-1))
allocate(g_k2(0:dim_hk-1,0:dim_hk-1))
allocate(g_k_res(0:dim_hk-1,0:dim_hk-1))

l_1l=.true.

do en=e_min,e_max,delta_e
dos_k=0
g_k=0
g_k2=0

do ik=0, nk_plot-1 	!k_point
g_k_res=0
    do i=0,dim_hk-1	!eigenvalue
!    weight=lorentz(en-e_k_tot_nscf(i,ik),T_dos)*k_vecs_nscf(ik,4)
!    weight=fermi(en-e_k_tot_nscf(i,ik),T_dos)*k_vecs_nscf(ik,4)
weight_cmplx=1.0d0/(en-e_k_tot_nscf(i,ik)+(0,1)*T_dos)/pi*k_vecs_nscf(ik,4)
weight=-aimag(weight_cmplx)

!dos
do j=0,dim_hk-1	!orbital
dos_k(j)=dos_k(j)+weight*abs(w_k_tot_nscf(j,i,ik))**2	!DOS
enddo
!23 DOS C 45 DOS F
!67 Dos cf im/re 
!89 dos cf kz>0 im/re

!df
do j=0,2*l_orbital-1,2	  !d	i just take down spin, otherwise it sums to zero
do k=2*l_orbital,dim_hk-1,2 !f
dos_k(dim_hk+0)=dos_k(dim_hk+0)-2*aimag(weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik)))		!df interference (2 for fd)
dos_k(dim_hk+1)=dos_k(dim_hk+1)+2*dble (weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik)))		!df interference (2 for fd)
if(k_vecs_nscf(ik,3)>=0) dos_k(dim_hk+2)=dos_k(dim_hk+2)-2*aimag(weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik)))		!df interference (2 for fd)
if(k_vecs_nscf(ik,3)>=0) dos_k(dim_hk+3)=dos_k(dim_hk+3)+2*dble (weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik)))	!df interference (2 for fd)
enddo
enddo

do j=0,dim_hk-1	!orbital
do k=0,dim_hk-1	!orbital
if(k_vecs_nscf(ik,3)>=0) g_k(j,k)=g_k(j,k)- weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik))
			g_k2(j,k)=g_k2(j,k)-weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik))
enddo
enddo

do j=0,dim_hk-1	!orbital
do k=0,dim_hk-1	!orbital
g_k_res(j,k)=g_k_res(j,k)-weight_cmplx*w_k_tot_nscf(j,i,ik)*conjg(w_k_tot_nscf(k,i,ik))/k_vecs_nscf(ik,4)
enddo
enddo


enddo	!i

!summed just over eigenvalue, not k
!if (l_1l) then
!write(700,"(i6,20f15.5)" ), ik, (k_vecs_nscf(ik,j),j=1,4)
!write(700,"(f10.5,20e15.5)" ), en
!do j=0,dim_hk-1	!orbital
!write(700,"(40f12.5)" ), (g_k_res(j,k),  k=0,dim_hk-1)	!orbital
!enddo
!write(700,"(a1)" ), ""
!	endif

enddo	!ik



write(600,"(f10.5,20e15.5)" ), en
do j=0,dim_hk-1	!orbital
write(600,"(40f12.5)" ), (g_k(j,k),  k=0,dim_hk-1)	!summed over k_z>0
enddo
write(600,"(a1)" ), ""

write(800,"(f10.5,20e15.5)" ), en
do j=0,dim_hk-1	!orbital
write(800,"(40f12.5)" ), (g_k2(j,k),  k=0,dim_hk-1)	!summed over all k
enddo
write(800,"(a1)" ), ""



write(500,"(f10.5,20e15.5)" ), en, (dos_k(j),j=0,dim_hk+3)

l_1l=.false.	!i do the k-resolved g just for the first energy for size reasons

enddo

 deallocate(dos_k)
 close(unit=500)
 close(unit=600)
endif

!    do i=0,size_h_kr-1
 !    do ik= 0, nkr_plot-1
  !    write(26,"(i5,3f10.5, 5f15.10)") ik, kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), kr_vecs_nscf(ik,3),e_kr_tot_nscf(i,ik)
   !  enddo
    ! write(26,"(a10)") "   "
    !enddo



!if (abs(T-T_min)<1e-6 .and. lspectrfunct) call spectralfunction_k_nscf(e_k_tot_nscf,w_k_tot_nscf)

!if (abs(T-T_min)<1e-6 .and. ldos) call compute_dos_k_nscf


!PATH 3d

write(*,"(a20)") "Nscf 3d"
print *, nk_plot_path3d

do ik=0, nk_plot_path3d-1

!!!builds and diag hamiltonian

  call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk_plot_path3d,phi_k_sigma_alpha_path3d,phi_k_alpha_sigma_path3d, kin_energy_k_path3d, k_vecs_path3d)

  call diagonalize_ham_k(Ham_k,e_k,w_k)

!writes ham and eigen
!  write (*, "(i5, 3f10.5)")ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2),k_vecs_path3d(ik,3)
!  call write_ham_k(ham_k)
!  do i=0,4*l_orbital-1
!  write (*, "(f15.7)") e_k(i)
!  enddo

!  e_k=e_k-mu_k
  e_k_tot_path3d(:,ik)=e_k(:)

ew_k_tot_path3d(:,ik,:)=0

do i=0,l_orb_tot-1
ew_k_tot_path3d(:,ik,i)=ew_k_tot_path3d(:,ik,i)+abs(w_k(2*i,:))**2+abs(w_k(2*i+1,:))**2	! weight on orbitals
enddo


enddo

write(24,"(a20)") "K points nscf 3d"
do ik=0, nk_plot_path3d-1 

 write(24,"(i6,6f10.5)")  ik, 	 k_vecs_path3d(ik,1), k_vecs_path3d(ik,2), k_vecs_path3d(ik,3),delta_k_path3d(ik)


 !    do is=0,1
!     write(24,"(4 ('(',f6.2,',',f6.2,')','  ') )")  real(phi_k_sigma_alpha_path3d(0,0,is,0,ik)),aimag(phi_k_sigma_alpha_path3d(0,0,is,0,ik)), &
 !    						   real(phi_k_sigma_alpha_path3d(0,0,is,1,ik)),aimag(phi_k_sigma_alpha_path3d(0,0,is,1,ik)),&
 !  						   real(phi_k_alpha_sigma_path3d(0,0,is,0,ik)),aimag(phi_k_alpha_sigma_path3d(0,0,is,0,ik)), &
 !  						   real(phi_k_alpha_sigma_path3d(0,0,is,1,ik)),aimag(phi_k_alpha_sigma_path3d(0,0,is,1,ik))
 !   enddo 

  write(24, " ")
enddo
write(24,"(a10)") "**********"

!bands
OPEN(unit=100,file=trim(label)//'/bands_k_3d',status='unknown')

  write(100,"(3a10, 5a15)") "#k_x","k_y","k_z", "E", "Weights"


do i=0,dim_hk-1
 do  ik=0, nk_plot_path3d-1
      write(100,"(i5,3f10.5, 10f15.10)") ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2), k_vecs_path3d(ik,3),e_k_tot_path3d(i,ik),(ew_k_tot_path3d(i,ik,j), j=0,l_orb_tot-1)
 enddo
  write(100,"(a10)") "   "
  write(100,"(a10)") "   "

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!k dot p
if (l_orbital==2 .or. phi_type=="sm2") then

write(*,"(a20)") "Nscf k dot p"
print *, nk_plot_path3d


allocate (ham_k_red(0:3,0:3))
allocate (e_k_red(0:3))
allocate (w_k_red(0:3,0:3))
allocate (e_k_tot_path3d_red(0:3,0:nk_plot_path3d-1))
allocate (rot_kdp(0:3,0:dim_hk-1))


  call set_param_kdotp
  call set_rot_kdotp(rot_kdp)

do ik=0, nk_plot_path3d-1

!!!builds and diag hamiltonian

  call build_ham_kdotp(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk_plot_path3d,phi_k_sigma_alpha_path3d,phi_k_alpha_sigma_path3d, kin_energy_k_path3d, k_vecs_path3d)
  call build_ham_kdotp_red(Ham_k,Ham_k_red,rot_kdp)

  CALL DIAGONALIZE( 'V', 4, 4, Ham_k_red, E_k_red, W_k_red )
  call diagonalize_ham_k(Ham_k,e_k,w_k)

!writes ham and eigen
!  write (*, "(i5, 3f10.5)")ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2),k_vecs_path3d(ik,3)
!  call write_ham_k(ham_k)
!  do i=0,4*l_orbital-1
!  write (*, "(f15.7)") e_k(i)
!  enddo

!  e_k=e_k-mu_k
  e_k_tot_path3d(:,ik)=e_k(:)
  e_k_tot_path3d_red(:,ik)=e_k_red(:)

ew_k_tot_path3d(:,ik,:)=0

do i=0,l_orb_tot-1
ew_k_tot_path3d(:,ik,i)=ew_k_tot_path3d(:,ik,i)+abs(w_k(2*i,:))**2+abs(w_k(2*i+1,:))**2	! weight on orbitals
enddo


enddo

write(24,"(a20)") "K points nscf 3d k dot p"
do ik=0, nk_plot_path3d-1 

 write(24,"(i6,6f10.5)")  ik, 	 k_vecs_path3d(ik,1), k_vecs_path3d(ik,2), k_vecs_path3d(ik,3),delta_k_path3d(ik)


  write(24, " ")
enddo
write(24,"(a10)") "**********"

!bands
OPEN(unit=100,file=trim(label)//'/bands_kdotp',status='unknown')

  write(100,"(3a10, 5a15)") "#k_x","k_y","k_z", "E", "Weights"


do i=0,dim_hk-1
 do  ik=0, nk_plot_path3d-1
      write(100,"(i5,3f10.5, 10f15.10)") ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2), k_vecs_path3d(ik,3),e_k_tot_path3d(i,ik),(ew_k_tot_path3d(i,ik,j), j=0,l_orb_tot-1)
 enddo
  write(100,"(a10)") "   "
  write(100,"(a10)") "   "

enddo

 close(unit=100)

OPEN(unit=100,file=trim(label)//'/bands_kdotp_red',status='unknown')

  write(100,"(3a10, 5a15)") "#k_x","k_y","k_z", "E", "Weights"


do i=0,3
 do  ik=0, nk_plot_path3d-1
      write(100,"(i5,3f10.5, 10f15.10)") ik, k_vecs_path3d(ik,1),k_vecs_path3d(ik,2), k_vecs_path3d(ik,3),e_k_tot_path3d_red(i,ik)
 enddo
  write(100,"(a10)") "   "
  write(100,"(a10)") "   "

enddo

 close(unit=100)
 
 
 deallocate (e_k_tot_path3d_red)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(a20)") "High symm points"

print *, nk_plot_hsp

print *, "0: Identity"
print *, "1: Parity"
print *, "2: TRS"
print *, "3: M_yz/i"
print *, "4: M_xz/i"
print *, "5: M_xy/i"
print *, "6: M_xz-yz/i"
print *, "7: M_xz+yz/i"
print *, "8: R_z pi/2"
print *, "9: R_z pi"
print *,"10: R_(x-y) pi"


 write(*,"(a14,30a13)"),"energy", "I", "iP", "TRS", "M_yz", "M_xz", "M_xy", "M_x-y", "M_x+y", "M_x-z", "M_x+z", "M_y-z", "M_y+z"


!  w_k_tot_hsp=0
!  e_k_tot_hsp=0

do ik=0, nk_plot_hsp-1

print *, ""
 write(*,"(a10,i6,6f10.5)") "**********",  ik, 	 k_vecs_hsp(ik,1), k_vecs_hsp(ik,2), k_vecs_hsp(ik,3)


!!!builds and diag hamiltonian

  call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk_plot_hsp,phi_k_sigma_alpha_hsp,phi_k_alpha_sigma_hsp, kin_energy_k_hsp, k_vecs_hsp)
  
  call diagonalize_ham_k(Ham_k,e_k,w_k)

 ! w_k_tot_hsp(:,ik,:)=w_k(:,:)
 ! e_k_tot_hsp(:,ik)=e_k(:)


if (ik==0) then
  n_occ=0	!number of occupied bands at gamma - I assume it is the same for all points
  n_empty=0	!number of empty bands at gamma
  do i=0,dim_hk-1
  if(e_k(i)<0) n_occ=n_occ+1
  enddo
  n_empty=dim_hk-n_occ

if (mod(n_occ,2)==1) print *, "warning, n_occ odd!"

allocate(parity_hsp_f(0:nk_plot_hsp-1,0:n_occ-1))
allocate(parity_hsp_e(0:nk_plot_hsp-1,0:n_empty-1))
parity_hsp_f=0
parity_hsp_e=0
endif

  print *, "n occ bands:", n_occ
  print *, "n empty bands:", n_empty

!on full states
print *, "filled states", n_occ

  allocate(symm_hsp(0:n_occ-1,0:n_occ-1))
  allocate(symm_hsp_p(0:n_occ/2-1,0:n_occ/2-1))
  allocate(symm_hsp_m(0:n_occ/2-1,0:n_occ/2-1))
  allocate(symm_hsp_pm(0:n_occ/2-1,0:n_occ/2-1))
  allocate(symm_hsp_mp(0:n_occ/2-1,0:n_occ/2-1))
  allocate(eigen_symm(0:n_occ/2-1))
  allocate(eigen_symm_real(0:n_occ-1))
  allocate(w_symm(0:n_occ-1,0:n_occ-1))
  allocate(w_symm_xy_pf(0:n_occ-1,0:n_occ/2-1))	!eigenvector for xy reflection  + full
  allocate(w_symm_xy_mf(0:n_occ-1,0:n_occ/2-1))	!eigenvector for xy reflection  - full

!symmetry matrix on occupied states
  do ixyz=0,n_sym_arpes
  
   symm_hsp=0
   eigen_symm=0
   eigen_symm_real=0

   do i=0, n_occ-1	!eigenvalue
   do j=0, n_occ-1	!eigenvalue
   do k=0, dim_hk-1
   do l=0, dim_hk-1

   symm_hsp(i,j)=symm_hsp(i,j)+conjg(w_k(k,i))*simm_oper(k,l,ixyz)*w_k(l,j)		!n_occ*n_occ!

   enddo
   enddo
   enddo
   enddo

   do i=0,n_occ-1
   write (*, "(30f10.5)"), (symm_hsp(i,j), j=0,n_occ-1)
   ENDDO



!(Anti)Hermitian operators
   if (ixyz<8) then
     call DIAGONALIZE( 'V', n_occ, n_occ, symm_hsp, eigen_symm_real, W_symm )!reflections: hermitian operator (times i, that I have removed)
     write(*,"(i3,f8.3, 30f7.2)") ixyz,(product(eigen_symm_real)), (eigen_symm_real(i), i=0,n_occ-1)

   if (ixyz==5) then	!M_xy reflection
    j=0
    do i=0, n_occ-1
     if    (eigen_symm_real(i)<0) then
      w_symm_xy_mf(:,j)=w_symm(:,i)
      j=j+1
     endif
    enddo
    j=0
    do i=0, n_occ-1
     if    (eigen_symm_real(i)>0) then
      w_symm_xy_pf(:,j)=w_symm(:,i)
      j=j+1
     endif
    enddo
   endif

   if (ixyz==1) parity_hsp_f(ik,:)=eigen_symm_real(:)

   else !rotations , either for states with + or - mirror eigenvalue
   

   symm_hsp_p = MATMUL( TRANSPOSE(real(w_symm_xy_pf)-(0,1)*aimag(w_symm_xy_pf)), MATMUL(symm_hsp,w_symm_xy_pf)  )
   symm_hsp_m = MATMUL( TRANSPOSE(real(w_symm_xy_mf)-(0,1)*aimag(w_symm_xy_mf)), MATMUL(symm_hsp,w_symm_xy_mf)  )
   symm_hsp_pm= MATMUL( TRANSPOSE(real(w_symm_xy_pf)-(0,1)*aimag(w_symm_xy_pf)), MATMUL(symm_hsp,w_symm_xy_mf)  )
   symm_hsp_mp= MATMUL( TRANSPOSE(real(w_symm_xy_mf)-(0,1)*aimag(w_symm_xy_mf)), MATMUL(symm_hsp,w_symm_xy_pf)  )

!print *, "pp"
!   do i=0,n_occ/2-1
!   write (*, "(30f10.5)"), (symm_hsp_p(i,j), j=0,n_occ/2-1)
!   ENDDO
!print *, "mm"
!   do i=0,n_occ/2-1
!   write (*, "(30f10.5)"), (symm_hsp_m(i,j), j=0,n_occ/2-1)
!   ENDDO
!print *, "pm"
!   do i=0,n_occ/2-1
!   write (*, "(30f10.5)"), (symm_hsp_pm(i,j), j=0,n_occ/2-1)
!   ENDDO
!print *, "mp"
!   do i=0,n_occ/2-1
!   write (*, "(30f10.5)"), (symm_hsp_mp(i,j), j=0,n_occ/2-1)
!   ENDDO

!   do i=0, n_occ/2-1
!   do j=0, n_occ/2-1
!    if(abs(symm_hsp_mp(i,j)>1e-5) print *, "WARNING, symm_hsp_mp not zero"
!    if(abs(symm_hsp_pm(i,j)>1e-5) print *, "WARNING, symm_hsp_mp not zero"
!   enddo
!   enddo
  
   call diagonalize_s(n_occ/2,symm_hsp_p,eigen_symm)	!rotations: non-hermitian operator
!   write(*,"(i3,a1,f8.3, f8.3,30(a2,f5.2,a1,f5.2,a1))") ixyz,"+",real(product(eigen_symm)), aimag(product(eigen_symm)), (" (",real(eigen_symm(i)),",",aimag(eigen_symm(i)),")", i=0,n_occ/2-1)

   call diagonalize_s(n_occ/2,symm_hsp_m,eigen_symm)	!rotations: non-hermitian operator
!   write(*,"(i3,a1,f8.3, f8.3,30(a2,f5.2,a1,f5.2,a1))") ixyz,"-",real(product(eigen_symm)), aimag(product(eigen_symm)), (" (",real(eigen_symm(i)),",",aimag(eigen_symm(i)),")", i=0,n_occ/2-1)

  endif

!plus and minus eigenvectors for reflection



  enddo	!ixyz
  deallocate(symm_hsp)
  deallocate(symm_hsp_p)
  deallocate(symm_hsp_m)
  deallocate(symm_hsp_pm)
  deallocate(symm_hsp_mp)
  deallocate(eigen_symm)
  deallocate(eigen_symm_real)
  deallocate(w_symm)


!empty states
 print *, "empty states", n_empty
  allocate(symm_hsp(0:n_empty-1,0:n_empty-1))
  allocate(symm_hsp_p(0:n_empty/2-1,0:n_empty/2-1))
  allocate(symm_hsp_m(0:n_empty/2-1,0:n_empty/2-1))
  allocate(eigen_symm(0:n_empty/2-1))
  allocate(eigen_symm_real(0:n_empty-1))
  allocate(w_symm(0:n_empty-1,0:n_empty-1))
  allocate(w_symm_xy_pe(0:n_empty-1,0:n_empty/2-1))	!eigenvector for xy reflection 1: - empty, 2: + empty
  allocate(w_symm_xy_me(0:n_empty-1,0:n_empty/2-1))	!eigenvector for xy reflection 1: - empty, 2: + empty

  do ixyz=0,n_sym_arpes
   symm_hsp=0
   eigen_symm=0
   eigen_symm_real=0

   do i=n_occ, dim_hk-1	!eigenvalue
   do j=n_occ, dim_hk-1	!eigenvalue
   do k=0, dim_hk-1
   do l=0, dim_hk-1

   symm_hsp(i-n_occ,j-n_occ)=symm_hsp(i-n_occ,j-n_occ)+conjg(w_k(k,i))*simm_oper(k,l,ixyz)*w_k(l,j)

   enddo
   enddo
   enddo
   enddo

   do i=0,n_empty-1
   write (*, "(30f10.5)"), (symm_hsp(i,j), j=0,n_empty-1)
   ENDDO



!(Anti)Hermitian operators
   if (ixyz<8) then
     call DIAGONALIZE( 'V', n_empty, n_empty, symm_hsp, eigen_symm_real, W_symm )!reflections: hermitian operator (times i, that I have removed)
     write(*,"(i3,f8.3, 30f7.2)") ixyz,(product(eigen_symm_real)), (eigen_symm_real(i), i=0,n_empty-1)

   if (ixyz==5) then	!M_xy reflection
    j=0
    do i=0, n_empty-1
     if    (eigen_symm_real(i)<0) then
      w_symm_xy_me(:,j)=w_symm(:,i)
      j=j+1
     endif
    enddo
    j=0
    do i=0, n_empty-1
     if    (eigen_symm_real(i)>0) then
      w_symm_xy_pe(:,j)=w_symm(:,i)
      j=j+1
     endif
    enddo
   endif

      if (ixyz==1) parity_hsp_e(ik,:)=eigen_symm_real(:)


   else !rotations , either for states with + or - mirror eigenvalue
   

   symm_hsp_p = MATMUL( TRANSPOSE(real(w_symm_xy_pe)-(0,1)*aimag(w_symm_xy_pe)), MATMUL(symm_hsp,w_symm_xy_pe)  )
   symm_hsp_m = MATMUL( TRANSPOSE(real(w_symm_xy_me)-(0,1)*aimag(w_symm_xy_me)), MATMUL(symm_hsp,w_symm_xy_me)  )

   
   call diagonalize_s(n_empty/2,symm_hsp_p,eigen_symm)	!rotations: non-hermitian operator
   write(*,"(i3,a1,f8.3, f8.3,30(a2,f5.2,a1,f5.2,a1))") ixyz,"+",real(product(eigen_symm)), aimag(product(eigen_symm)), (" (",real(eigen_symm(i)),",",aimag(eigen_symm(i)),")", i=0,n_empty/2-1)

   call diagonalize_s(n_empty/2,symm_hsp_m,eigen_symm)	!rotations: non-hermitian operator
   write(*,"(i3,a1,f8.3, f8.3,30(a2,f5.2,a1,f5.2,a1))") ixyz,"-",real(product(eigen_symm)), aimag(product(eigen_symm)), (" (",real(eigen_symm(i)),",",aimag(eigen_symm(i)),")", i=0,n_empty/2-1)

  endif

!plus and minus eigenvectors for reflection



  enddo	!ixyz
  
  
  
  
  
  
  deallocate(symm_hsp)
  deallocate(symm_hsp_p)
  deallocate(symm_hsp_m)
  deallocate(eigen_symm)
  deallocate(eigen_symm_real)
  deallocate(w_symm)
  deallocate(w_symm_xy_pe)
  deallocate(w_symm_xy_me)
  deallocate(w_symm_xy_pf)
  deallocate(w_symm_xy_mf)


 enddo	!ik

!do i=0, dim_hk-1
! write(*,"(i1,f13.8,30(a1,f5.2,a1,f5.2,a1))") i,e_k(i), (("(",real(symm_hsp(ixyz,i)),",",aimag(symm_hsp(ixyz,i)),")"),ixyz=0,n_sym_arpes)
!enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Chern number
if (lchern) then
print *, "*************************"
print *, "Chern number for M_yz/xz/xy and M_xz-yz"

allocate(M_xy(0:dim_hk-1,0:dim_hk-1))
if (iplane_xyz==1) M_xy(:,:)=simm_oper(:,:,3)		!3:M_yz	4:M_xz 5:M_xy
if (iplane_xyz==2) M_xy(:,:)=simm_oper(:,:,4)		!3:M_yz	4:M_xz 5:M_xy
if (iplane_xyz==3) M_xy(:,:)=simm_oper(:,:,5)		!3:M_yz	4:M_xz 5:M_xy
if (iplane_xyz==1) print *, "Reflection symmetry operator M_yz on basis"
if (iplane_xyz==2) print *, "Reflection symmetry operator M_xz on basis"
if (iplane_xyz==3) print *, "Reflection symmetry operator M_xy on basis"
do i=0,dim_hk-1
 write(*,"(30(a1,f5.2,a1,f5.2,a2))") ( "(",real(m_xy(i,j)),",",aimag(m_xy(i,j)),") ", j=0,dim_hk-1 )
enddo


allocate(e_xy(0:dim_hk-1))
allocate(w_xy(0:dim_hk-1,0:dim_hk-1))
allocate(w_xy_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_xy_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(e_xzyz(0:dim_hk-1))
allocate(w_xzyz(0:dim_hk-1,0:dim_hk-1))
allocate(w_xzyz_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_xzyz_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(M_check(0:dim_hk-1,0:dim_hk-1))
allocate(M_check2(0:dim_hk/2-1,0:dim_hk/2-1))


!M_xy
call DIAGONALIZE( 'V', dim_hk, dim_hk, M_xy, e_xy, w_xy)!reflections: hermitian operator (times i, that I have removed)

do i=0, dim_hk/2-1
w_xy_m(:,i)=w_xy(:,i)		!eigenstates -
w_xy_p(:,i)=w_xy(:,i+dim_hk/2)	!eigenstates +
enddo


print *, "eigenvectors of M_yz/xz/xy:"
    write(*,"(100f15.9)") ( e_xy(j), j=0,dim_hk-1 )
do, i=0,dim_hk-1
    write(*,"(100(a1,f6.2,a1,f6.2,a1))") ( "(",real(w_xy(i,j)),",",aimag(w_xy(i,j)),")", j=0,dim_hk-1 )
enddo


m_check=0
call basis_rotationcmplx_c(w_xy,M_xy,w_xy, M_check,dim_hk,dim_hk)
print *, "M_xy diag"
do i=0,dim_hk-1
 write(*,"(30(a1,f5.2,a1,f5.2,a2))") ( "(",real(m_check(i,j)),",",aimag(m_check(i,j)),") ", j=0,dim_hk-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xy_p,M_xy,w_xy_p, m_check2,dim_hk,dim_hk/2)
print *, "M++"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xy_m,M_xy,w_xy_m, m_check2,dim_hk,dim_hk/2)
print *, "M--"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xy_p,M_xy,w_xy_m, m_check2,dim_hk,dim_hk/2)
print *, "M+-"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xy_m,M_xy,w_xy_p, m_check2,dim_hk,dim_hk/2)
print *, "M-+"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

!M_xzyz
allocate(M_xzyz(0:dim_hk-1,0:dim_hk-1))
M_xzyz(:,:)=simm_oper(:,:,6)
print *, "Reflection symmetry operator M_xz-yz on basis"
do i=0,dim_hk-1
 write(*,"(30(a1,f5.2,a1,f5.2,a2))") ( "(",real(m_xzyz(i,j)),",",aimag(m_xzyz(i,j)),") ", j=0,dim_hk-1 )
enddo

call DIAGONALIZE( 'V', dim_hk, dim_hk, M_xzyz, e_xzyz, w_xzyz)!reflections: hermitian operator (times i, that I have removed)

do i=0, dim_hk/2-1
w_xzyz_m(:,i)=w_xzyz(:,i)		!eigenstates -
w_xzyz_p(:,i)=w_xzyz(:,i+dim_hk/2)	!eigenstates +
enddo

print *, "eigenvectors of M_xzyz:"
    write(*,"(100f15.9)") ( e_xzyz(j), j=0,dim_hk-1 )
do, i=0,dim_hk-1
    write(*,"(100(a1,f6.2,a1,f6.2,a1))") ( "(",real(w_xzyz(i,j)),",",aimag(w_xzyz(i,j)),")", j=0,dim_hk-1 )
enddo


m_check=0
call basis_rotationcmplx_c(w_xzyz,M_xzyz,w_xzyz, m_check,dim_hk,dim_hk)
print *, "M_xzyz diag"
do i=0,dim_hk-1
 write(*,"(30(a1,f5.2,a1,f5.2,a2))") ( "(",real(m_check(i,j)),",",aimag(m_check(i,j)),") ", j=0,dim_hk-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xzyz_p,M_xzyz,w_xzyz_p, m_check2,dim_hk,dim_hk/2)
print *, "M++"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xzyz_m,M_xzyz,w_xzyz_m, m_check2,dim_hk,dim_hk/2)
print *, "M--"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xzyz_p,M_xzyz,w_xzyz_m, m_check2,dim_hk,dim_hk/2)
print *, "M+-"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo

m_check2=0
call basis_rotationcmplx_c(w_xzyz_m,M_xzyz,w_xzyz_p, m_check2,dim_hk,dim_hk/2)
print *, "M-+"
do i=0,dim_hk/2-1
 write(*,"(30(a1,f10.7,a1,f10.7,a2))") ( "(",real(m_check2(i,j)),",",aimag(m_check2(i,j)),") ", j=0,dim_hk/2-1 )
enddo


deallocate(w_xy)
deallocate(e_xy)
deallocate(w_xzyz)
deallocate(e_xzyz)

allocate(det_(4))

allocate(ham_k_p(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(ham_k_m(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(ham_k_pm(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(ham_k_mp(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(ham_k_tmp(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_k_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_k_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_k_p_tmp(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(w_k_m_tmp(0:dim_hk/2-1,0:dim_hk/2-1))
allocate(e_k_p(0:dim_hk/2-1))
allocate(e_k_m(0:dim_hk/2-1))
allocate(w_m_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_m_m(0:dim_hk-1,0:dim_hk/2-1))


do iplane=1,3

call set_k_vecs_chern(iplane)
if(iplane==1 .and. iplane_xyz==1) print *, "k_x=0"
if(iplane==1 .and. iplane_xyz==2) print *, "k_y=0"
if(iplane==1 .and. iplane_xyz==3) print *, "k_z=0"
if(iplane==2 .and. iplane_xyz==1) print *, "k_x=pi"
if(iplane==2 .and. iplane_xyz==2) print *, "k_y=pi"
if(iplane==2 .and. iplane_xyz==3) print *, "k_z=pi"
if(iplane==3) print *, "k_x=k_y"


allocate(e_k_tot_chern_p(0:dim_hk/2-1,0:nk_chern-1))
allocate(e_k_tot_chern_m(0:dim_hk/2-1,0:nk_chern-1))




if(iplane==1 .or. iplane==2) w_m_p=w_xy_p	!reflection in xy plane for kz=0 and kz=pi
if(iplane==3) w_m_p=w_xzyz_p			!reflection in xz-yz plane for kx=ky
if(iplane==1 .or. iplane==2) w_m_m=w_xy_m
if(iplane==3) w_m_m=w_xzyz_m


do ik=0, nk_chern-1
 !write(*, "(i8,30f10.5)"), ik,k_vecs_chern(ik,1), k_vecs_chern(ik,2), k_vecs_chern(ik,3)

  call build_ham_k(ik,Ham_k,b_k,lambda_k, mu_k,.true.,nk_chern,phi_k_sigma_alpha_chern,phi_k_alpha_sigma_chern, kin_energy_k_chern, k_vecs_chern)

!H++
ham_k_tmp=0
ham_k_p=0
call matmul2 (dim_hk,dim_hk/2,dim_hk, ham_k,w_m_p, ham_k_tmp) 	!M |psi>	!matmul2(M,N,K, A, B, X)  A( 0:M-1, 0:K-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
call matmul2hc2 (dim_hk/2,dim_hk/2,dim_hk, w_m_p, ham_k_tmp, ham_k_p) 	!<psi|(M |psi>) !matmul2hc2(M,N,K, A, B, X) A( 0:K-1, 0:M-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
!call basis_rotationcmplx_c(w_m_p,ham_k,w_m_p, ham_k_p,dim_hk,dim_hk/2)
!H--
ham_k_tmp=0
ham_k_m=0
call matmul2 (dim_hk,dim_hk/2,dim_hk, ham_k, w_m_m, ham_k_tmp) 	!M |psi>	!matmul2(M,N,K, A, B, X)  A( 0:M-1, 0:K-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
call matmul2hc2 (dim_hk/2,dim_hk/2,dim_hk, w_m_m, ham_k_tmp, ham_k_m) 	!<psi|(M |psi>) !matmul2hc2(M,N,K, A, B, X) A( 0:K-1, 0:M-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
!call basis_rotationcmplx_c(w_m_m,ham_k,w_m_m, ham_k_m,dim_hk,dim_hk/2)
!H+-
ham_k_tmp=0
ham_k_pm=0
call matmul2 (dim_hk,dim_hk/2,dim_hk, ham_k, w_m_m, ham_k_tmp) 	!M |psi>	!matmul2(M,N,K, A, B, X)  A( 0:M-1, 0:K-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
call matmul2hc2 (dim_hk/2,dim_hk/2,dim_hk, w_m_p, ham_k_tmp, ham_k_pm) 	!<psi|(M |psi>) !matmul2hc2(M,N,K, A, B, X) A( 0:K-1, 0:M-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
!call basis_rotationcmplx_c(w_m_p,ham_k,w_m_m, ham_k_pm,dim_hk,dim_hk/2)
!H-+
ham_k_tmp=0
ham_k_mp=0
call matmul2 (dim_hk,dim_hk/2,dim_hk, ham_k, w_m_p, ham_k_tmp) 	!M |psi>	!matmul2(M,N,K, A, B, X)  A( 0:M-1, 0:K-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
call matmul2hc2 (dim_hk/2,dim_hk/2,dim_hk, w_m_m, ham_k_tmp, ham_k_mp) 	!<psi|(M |psi>) !matmul2hc2(M,N,K, A, B, X) A( 0:K-1, 0:M-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
!call basis_rotationcmplx_c(w_m_m,ham_k,w_m_p, ham_k_mp,dim_hk,dim_hk/2)

!checks that H+- and H-+ are zero
do i=0, dim_hk/2-1
do j=0, dim_hk/2-1
!if (abs(ham_k_pm(i,j))>1e-4) write(*, "(a15,i10,3i5,2f15.10)") "WARNING H+-",ik,iplane,i,j,ham_k_pm(i,j)
!if (abs(ham_k_mp(i,j))>1e-4) write(*, "(a15,i10,3i5,2f15.10)") "WARNING H-+",ik,iplane,i,j,ham_k_mp(i,j)
enddo
enddo


!diag H++ and H--
call DIAGONALIZE( 'V', dim_hk/2, dim_hk/2, ham_k_p, e_k_p, w_k_p_tmp)	!H++
call DIAGONALIZE( 'V', dim_hk/2, dim_hk/2, ham_k_m, e_k_m, w_k_m_tmp)	!H--


!if (iplane==3 .and. abs(k_vecs_chern(ik,1))<1e-5) then
!write(*, "(i8,30f10.5)"), ik,k_vecs_chern(ik,1), k_vecs_chern(ik,2), k_vecs_chern(ik,3)
!do, i=0,dim_hk/2-1
!    write(*,"(100(a1,f15.10,a1,f15.10,a1))") ( "(",real(Ham_k(i,j)),",",aimag(Ham_k(i,j)),")", j=0,dim_hk/2-1 )
!enddo
!endif


if (ik==0) then
allocate(w_k_tot_chern_fp(0:dim_hk/2-1, 0:n_occ/2-1,0:nk_chern-1))
allocate(w_k_tot_chern_fm(0:dim_hk/2-1, 0:n_occ/2-1,0:nk_chern-1))
allocate(w_k_tot_chern_ep(0:dim_hk/2-1, 0:n_empty/2-1,0:nk_chern-1))
allocate(w_k_tot_chern_em(0:dim_hk/2-1, 0:n_empty/2-1,0:nk_chern-1))
endif

e_k_tot_chern_p(:,ik)=e_k_p(:)
e_k_tot_chern_m(:,ik)=e_k_m(:)

do i=0,n_occ/2-1
w_k_tot_chern_fm(:,i,ik)=w_k_m_tmp(:,i)
w_k_tot_chern_fp(:,i,ik)=w_k_p_tmp(:,i)
enddo
do i=0,n_empty/2-1
w_k_tot_chern_em(:,i,ik)=w_k_m_tmp(:,i+n_occ/2)
w_k_tot_chern_ep(:,i,ik)=w_k_p_tmp(:,i+n_occ/2)
enddo

enddo	!ik


!!!!!!!!U field
!now I have eigenvectors of Hamiltonian and M_xy, and I calculate the U field for + and - eigenstates of M_xy, and for filled and empty states
allocate(u_k_chern(0:nk_chern-1,2,4))	!field
allocate(u_mat_fm(0:n_occ/2-1,0:n_occ/2-1))
allocate(u_mat_fp(0:n_occ/2-1,0:n_occ/2-1))
allocate(u_mat_em(0:n_empty/2-1,0:n_empty/2-1))
allocate(u_mat_ep(0:n_empty/2-1,0:n_empty/2-1))

do ik=0, nk_chern-1
if (iplane==1 .or. iplane==2) call ik_pxpy (ik,nkx_chern, nky_chern, ik_px, ik_py)	!given ik, it gives ik_px just on the right, and ik_py just above
if (iplane==3) call ik_pxpy (ik,nkxy_chern, nkz_chern, ik_px, ik_py)	!given ik, it gives ik_px just on the right, and ik_py just above

do i=0,n_occ/2-1
 norm=0
 do k=0, dim_hk/2-1
  norm=norm+abs(w_k_tot_chern_fm(k,i,ik))**2
 enddo
 if (abs(norm-1)>1e-5) print *,"norm fm", ik, norm
 norm=0
 do k=0, dim_hk/2-1
  norm=norm+abs(w_k_tot_chern_fp(k,i,ik))**2
 enddo
 if (abs(norm-1)>1e-5) print *,"norm fp", ik, norm
enddo
do i=0,n_empty/2-1
 norm=0
 do k=0, dim_hk/2-1
  norm=norm+abs(w_k_tot_chern_em(k,i,ik))**2
 enddo
 if (abs(norm-1)>1e-5) print *,"norm em", ik, norm
 norm=0
 do k=0, dim_hk/2-1
  norm=norm+abs(w_k_tot_chern_ep(k,i,ik))**2
 enddo
 if (abs(norm-1)>1e-5) print *,"norm ep", ik, norm
enddo


!Field along x
u_mat_em=0
u_mat_ep=0
u_mat_fm=0
u_mat_fp=0

do k=0, dim_hk/2-1	!in the space of + states
do i=0,n_occ/2-1	!1/2 of them is filled
do j=0,n_occ/2-1
u_mat_fm(i,j)=u_mat_fm(i,j)+conjg(w_k_tot_chern_fm(k,i,ik))*w_k_tot_chern_fm(k,j,ik_px)
u_mat_fp(i,j)=u_mat_fp(i,j)+conjg(w_k_tot_chern_fp(k,i,ik))*w_k_tot_chern_fp(k,j,ik_px)
enddo
enddo
enddo

do k=0, dim_hk/2-1
do i=0,n_empty/2-1
do j=0,n_empty/2-1
u_mat_em(i,j)=u_mat_em(i,j)+conjg(w_k_tot_chern_em(k,i,ik))*w_k_tot_chern_em(k,j,ik_px)
u_mat_ep(i,j)=u_mat_ep(i,j)+conjg(w_k_tot_chern_ep(k,i,ik))*w_k_tot_chern_ep(k,j,ik_px)
enddo
enddo
enddo


!with det
det_(1)=det_man(n_occ/2, u_mat_fm)
det_(2)=det_man(n_occ/2, u_mat_fp)
det_(3)=det_man(n_empty/2, u_mat_em)
det_(4)=det_man(n_empty/2, u_mat_ep)

u_k_chern(ik,1,:)=det_(:)/abs(det_(:))

if(abs(abs(det_(1))-1)>0.2) then
!do i=0,n_occ/2-1
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(u_mat_fm(i,j)),",",aimag(u_mat_fm(i,j)),") ", j=0,n_occ/2-1 )
!enddo
write(*, "(a3,i8,i8,6f10.5,10f15.8)"), "x",ik,ik_px, k_vecs_chern(ik,1),k_vecs_chern(ik,2),k_vecs_chern(ik,3),&
					      k_vecs_chern(ik_px,1),k_vecs_chern(ik_px,2),k_vecs_chern(ik_px,3),abs(det_(1))
!do i=0,n_occ/2-1
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(w_k_tot_chern_fm(k,i,ik)),",",aimag(w_k_tot_chern_fm(k,i,ik)),") ", k=0,dim_hk/2-1 )		
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(w_k_tot_chern_fm(k,i,ik_px)),",",aimag(w_k_tot_chern_fm(k,i,ik_px)),") ", k=0,dim_hk/2-1 )
! enddo
endif


!

!Field along y
u_mat_em=0
u_mat_ep=0
u_mat_fm=0
u_mat_fp=0


do k=0, dim_hk/2-1
do i=0,n_occ/2-1
do j=0,n_occ/2-1
u_mat_fm(i,j)=u_mat_fm(i,j)+conjg(w_k_tot_chern_fm(k,i,ik))*w_k_tot_chern_fm(k,j,ik_py)
u_mat_fp(i,j)=u_mat_fp(i,j)+conjg(w_k_tot_chern_fp(k,i,ik))*w_k_tot_chern_fp(k,j,ik_py)
enddo
enddo
enddo

do k=0, dim_hk/2-1
do i=0,n_empty/2-1
do j=0,n_empty/2-1
u_mat_em(i,j)=u_mat_em(i,j)+conjg(w_k_tot_chern_em(k,i,ik))*w_k_tot_chern_em(k,j,ik_py)
u_mat_ep(i,j)=u_mat_ep(i,j)+conjg(w_k_tot_chern_ep(k,i,ik))*w_k_tot_chern_ep(k,j,ik_py)
enddo
enddo
enddo



!with det
det_(1)=det_man(n_occ/2, u_mat_fm)
det_(2)=det_man(n_occ/2, u_mat_fp)
det_(3)=det_man(n_empty/2, u_mat_em)
det_(4)=det_man(n_empty/2, u_mat_ep)

u_k_chern(ik,2,:)=det_(:)/abs(det_(:))

if(abs(abs(det_(1))-1)>0.2) then
!do i=0,n_occ/2-1
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(u_mat_fm(i,j)),",",aimag(u_mat_fm(i,j)),") ", j=0,n_occ/2-1 )
!enddo
write(*, "(a3,i8,i8,6f10.5,10f15.8)"), "y",ik,ik_px, k_vecs_chern(ik,1),k_vecs_chern(ik,2),k_vecs_chern(ik,3),&
					      k_vecs_chern(ik_py,1),k_vecs_chern(ik_py,2),k_vecs_chern(ik_py,3),abs(det_(1))
!do i=0,n_occ/2-1
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(w_k_tot_chern_fm(k,i,ik)),",",aimag(w_k_tot_chern_fm(k,i,ik)),") ", k=0,dim_hk/2-1 )		
! write(*,"(30(a1,f15.10,a1,f15.10,a2))") ( "(",real(w_k_tot_chern_fm(k,i,ik_py)),",",aimag(w_k_tot_chern_fm(k,i,ik_py)),") ", k=0,dim_hk/2-1 )
! enddo
endif



!write(*, "(i8,10f15.8)"), ik,atan2(aimag(u_k_chern(ik,1,1)), real(u_k_chern(ik,1,1))),atan2(aimag(u_k_chern(ik,2,1)), real(u_k_chern(ik,2,1)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!print *, "u_k_x, u_k_y"
!write(*, "(20f10.5)") ,u_k_chern(ik,1,1), u_k_chern(ik,2,1)
!write(*, "(20f10.5)") ,u_k_chern(ik,1,2), u_k_chern(ik,2,2)
!write(*, "(20f10.5)") ,u_k_chern(ik,1,3), u_k_chern(ik,2,3)
!write(*, "(20f10.5)") ,u_k_chern(ik,1,4), u_k_chern(ik,2,4)
enddo	!ik





deallocate(u_mat_fp)
deallocate(u_mat_ep)
deallocate(u_mat_fm)
deallocate(u_mat_em)
!deallocate(e_u_fp)
!deallocate(e_u_ep)
!deallocate(e_u_fm)
!deallocate(e_u_em)
deallocate(w_k_tot_chern_fp)
deallocate(w_k_tot_chern_fm)
deallocate(w_k_tot_chern_ep)
deallocate(w_k_tot_chern_em)

print *, "f"
!Now I compute the Field F12 ...
allocate(f12_k(0:nk_chern-1,1:4))

do ik=0, nk_chern-1
if (iplane==1 .or. iplane==2) call ik_pxpy (ik,nkx_chern, nky_chern, ik_px, ik_py)	!given ik, it gives ik_px just on the right, and ik_py just above
if (iplane==3) call ik_pxpy (ik,nkxy_chern, nkz_chern, ik_px, ik_py)	!given ik, it gives ik_px just on the right, and ik_py just above

f12_k(ik,:)=log(u_k_chern(ik,1,:)*u_k_chern(ik_px,2,:)/u_k_chern(ik_py,1,:)/u_k_chern(ik,2,:))

enddo


!writes berry phase on file
if (iplane==1) OPEN(unit=100,file=trim(label)//'/berry_z0',status='unknown')
if (iplane==2) OPEN(unit=100,file=trim(label)//'/berry_zp',status='unknown')
if (iplane==3) OPEN(unit=100,file=trim(label)//'/berry_xy',status='unknown')
do ik=0, nk_chern-1
write (100, "(3f10.5,20e17.8)") k_vecs_chern(ik,1),k_vecs_chern(ik,2),k_vecs_chern(ik,3), (f12_k(ik,j),j=1,4)
enddo
 close(unit=100)


!... and the chern number
Write(*, "(a20,f10.5)"), "***********"
if(iplane==1 .and. iplane_xyz==1) Write(*, "(a20,f10.5)"), "k_x=0"
if(iplane==1 .and. iplane_xyz==2) Write(*, "(a20,f10.5)"), "k_y=0"
if(iplane==1 .and. iplane_xyz==3) Write(*, "(a20,f10.5)"), "k_z=0"
if(iplane==2 .and. iplane_xyz==1) Write(*, "(a20,f10.5)"), "k_x=pi"
if(iplane==2 .and. iplane_xyz==2) Write(*, "(a20,f10.5)"), "k_y=pi"
if(iplane==2 .and. iplane_xyz==3) Write(*, "(a20,f10.5)"), "k_z=pi"
Write(*, "(a20,f10.5)"), "Chern number f-:", -aimag(sum(f12_k(:,1))/(2*pi))	!-=i/2pi	
Write(*, "(a20,f10.5)"), "Chern number f+:", -aimag(sum(f12_k(:,2))/(2*pi))
Write(*, "(a20,f10.5)"), "Chern number e-:", -aimag(sum(f12_k(:,3))/(2*pi))
Write(*, "(a20,f10.5)"), "Chern number e+:", -aimag(sum(f12_k(:,4))/(2*pi))
Write(*, "(a20,f10.5)"), "***********"


OPEN(unit=300,file=trim(label)//'/indices',status='unknown')
if(iplane==1 .and. iplane_xyz==1) Write(300, "(a20,f10.5)"), "k_x=0"
if(iplane==1 .and. iplane_xyz==2) Write(300, "(a20,f10.5)"), "k_y=0"
if(iplane==1 .and. iplane_xyz==3) Write(300, "(a20,f10.5)"), "k_z=0"
if(iplane==2 .and. iplane_xyz==1) Write(300, "(a20,f10.5)"), "k_x=pi"
if(iplane==2 .and. iplane_xyz==2) Write(300, "(a20,f10.5)"), "k_y=pi"
if(iplane==2 .and. iplane_xyz==3) Write(300, "(a20,f10.5)"), "k_z=pi"
if(iplane==3) Write(300, "(a20,f10.5)"), "k_x=k_y"
Write(300, "(a20,f10.5)"), "Chern number f-:", -aimag(sum(f12_k(:,1))/(2*pi))
Write(300, "(a20,f10.5)"), "Chern number f+:", -aimag(sum(f12_k(:,2))/(2*pi))
Write(300, "(a20,f10.5)"), "Chern number e-:", -aimag(sum(f12_k(:,3))/(2*pi))
Write(300, "(a20,f10.5)"), "Chern number e+:", -aimag(sum(f12_k(:,4))/(2*pi))


!write energies for check of gap closing
if (iplane==1) OPEN(unit=1028,file=trim(label)//'/bands_k_chern_m1',status='unknown')
if (iplane==1) OPEN(unit=1029,file=trim(label)//'/bands_k_chern_p1',status='unknown')
if (iplane==2) OPEN(unit=1028,file=trim(label)//'/bands_k_chern_m2',status='unknown')
if (iplane==2) OPEN(unit=1029,file=trim(label)//'/bands_k_chern_p2',status='unknown')
if (iplane==3) OPEN(unit=1028,file=trim(label)//'/bands_k_chern_m3',status='unknown')
if (iplane==3) OPEN(unit=1029,file=trim(label)//'/bands_k_chern_p3',status='unknown')

 write(28,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,dim_hk/2-1
     do ik= 0, nk_chern-1
write(1028,"(i5,3f10.5, f17.10,15e12.5)") ik, k_vecs_chern(ik,1),k_vecs_chern(ik,2), k_vecs_chern(ik,3),e_k_tot_chern_m(i,ik)
write(1029,"(i5,3f10.5, f17.10,15e12.5)") ik, k_vecs_chern(ik,1),k_vecs_chern(ik,2), k_vecs_chern(ik,3),e_k_tot_chern_p(i,ik)
     enddo
write(1028,"(a10)") "   "
write(1028,"(a10)") "   "
write(1029,"(a10)") "   "
write(1029,"(a10)") "   "
    enddo
 close(unit=1028)
 close(unit=1029)



deallocate(e_k_tot_chern_m)
deallocate(e_k_tot_chern_p)


deallocate(f12_k)
deallocate(u_k_chern)
deallocate(k_vecs_chern)
deallocate(phi_k_sigma_alpha_chern)
deallocate(phi_k_alpha_sigma_chern)
deallocate(kin_energy_k_chern)
deallocate(delta_k_chern)
!deallocate(M_xy)
!deallocate(M_xy_e)
!deallocate(w_k_f)
!deallocate(w_k_e)
!deallocate(e_xy_f)
!deallocate(e_xy_e)
!deallocate(w_xy_f_tmp)
!deallocate(w_xy_e_tmp)
!deallocate(w_xy_f)
!deallocate(w_xy_e)
!deallocate(M_xy_tmp_f)
!deallocate(M_xy_tmp_e)
enddo	!iplane


!and I put topological indices, too
!  write (*, "(a40)") "Topological indices from parity:"

!print *, "Parities full:"
!do ik=0,7
!write(*, "(i5,20f10.5)"),ik,(parity_hsp_f(ik,j),j=0,n_occ-1)
!enddo
!print *, "Parities empty:"
!do ik=0,7
!write(*, "(i5,20f10.5)"),ik,(parity_hsp_e(ik,j),j=0,n_empty-1)
!enddo

!print *, "Indices full:" 
!call index_01(parity_hsp_f,n_occ,0)
!call index_01(parity_hsp_f,n_occ,1)
!call index_01(parity_hsp_f,n_occ,2)
!call index_01(parity_hsp_f,n_occ,3)

!print *, "Indices empty:"
!call index_01(parity_hsp_e,n_empty,0)
!call index_01(parity_hsp_e,n_empty,1)
!call index_01(parity_hsp_e,n_empty,2)
!call index_01(parity_hsp_e,n_empty,3)

 close(unit=300)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(k_vecs_nscf)
!deallocate(factor_k_nscf)
if (allocated(phi_k_sigma_alpha_nscf))deallocate(phi_k_sigma_alpha_nscf)
if (allocated(phi_k_alpha_sigma_nscf))deallocate(phi_k_alpha_sigma_nscf)
if (allocated(kin_energy_k_nscf))deallocate(kin_energy_k_nscf)
if (allocated(delta_k_nscf))deallocate(delta_k_nscf)
if (allocated(k_vecs_path3d))deallocate(k_vecs_path3d)
!deallocate(factor_k_path3d)
if (allocated(phi_k_sigma_alpha_path3d))deallocate(phi_k_sigma_alpha_path3d)
if (allocated(phi_k_alpha_sigma_path3d))deallocate(phi_k_alpha_sigma_path3d)
if (allocated(kin_energy_k_path3d))deallocate(kin_energy_k_path3d)
if (allocated(delta_k_path3d))deallocate(delta_k_path3d)
!deallocate(corr_fct_k_nscf)
deallocate(V_k_nscf)
!deallocate(e_corr_fct_k_nscf)
deallocate(fermien_k)
deallocate(occup_k)
if (allocated(ew_k_tot_nscf))deallocate(ew_k_tot_nscf)
if (allocated(w_k_tot_nscf))deallocate(w_k_tot_nscf)
if (allocated(ew_k_tot_path3d))deallocate(ew_k_tot_path3d)
deallocate(k_vecs_hsp)
!deallocate(factor_k_path3d)
deallocate(phi_k_sigma_alpha_hsp)
deallocate(phi_k_alpha_sigma_hsp)
deallocate(kin_energy_k_hsp)
deallocate(delta_k_hsp)
deallocate(M_xy)
deallocate(w_xy_p)
deallocate(w_xy_m)
deallocate(w_xzyz_p)
deallocate(w_xzyz_m)
deallocate(M_check)
deallocate(M_check2)

endif !lchern

end subroutine do_k_nscf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kvecs_kr_pw (k_vecs2d, nk2d,nkx, nky, icase, lkshifted)
real(kind=idp),allocatable, intent(out)	::k_vecs2d(:,:)
integer, intent(out)	::nk2d
integer, intent(in)	::nkx, nky
integer	::ik,ikn
integer,intent(in)	::icase
real	::dummy
 character(50)	::file_input,file_info, file_output
real(kind=idp)	:: y_corr, z_corr
logical, intent(in)	::lkshifted
inTEGER	::KSHIFTED




if (lysameasx .and. nx==ny) then
y_corr=1.0d0
else
y_corr=0.9999d0
endif
if (lzsameasx) then
z_corr=1.0d0
else
z_corr=0.9998d0
endif


if (lkshifted) kshifted=1
if (.not. lkshifted) kshifted=0

!0:kr
!1:r
if (icase==0) then
file_input="input_k2d"
file_info="info_k2d"
file_output="output_k2d"
else if (icase==1) then
file_input="input_r_k"
file_info="info_r_k"
file_output="output_r_k"
else if (icase==3) then
file_input="input_r_nscf"
file_info="info_r_nscf"
file_output="output_r_nscf"
endif

OPEN(unit=100,file=trim(label)//'/'//trim(file_input),status='unknown')


write(100,"(i1)") 8
write(100,"(a)") "k2d"
write(100,"(f15.8)") real(nkz)/real(nkx)*z_corr	!to impose that z direction is always different from x-y
write(100,"(f15.8)") real(nky)/real(nkx)*y_corr
write(100,"(3i6)") nkx,nky,1
write(100,"(3i6)") kshifted,kshifted,0
write(100,"(a1)") "f"

 CLOSE(unit=100)



 call system("kpoints.x <" // trim(label) // "/"// trim(file_input))
 call system("mv info " // trim(label)// "/" //trim(file_info) )
 call system("mv k2d " // trim(label)// "/" //trim(file_output) )


OPEN(unit=101,file=trim(label)//"/"// trim(file_output),status='unknown')
read(101,"(i)"), nk2d

allocate(k_vecs2d(0:nk2d-1,4))

do ik=0, nk2d-1
!read(101,"(i5,1x,3f20.10,f7.2)"), ikn, k_vecs2d(ikn-1,1),k_vecs2d(ikn-1,2),k_vecs2d(ikn-1,3),k_vecs2d(ikn-1,4)
read(101,*), ikn, k_vecs2d(ikn-1,1),k_vecs2d(ikn-1,2),k_vecs2d(ikn-1,3),k_vecs2d(ikn-1,4)
enddo
 CLOSE(unit=101)

write(*,"(a21,i5)") "Nk points from kpoints.x:", nk2d

do ik=0, nk2d-1
k_vecs2d(ik,2)=k_vecs2d(ik,2)*real(nky)/real(nkx)*y_corr
k_vecs2d(ik,3)=k_vecs2d(ik,3)*real(nkz)/real(nkx)*z_corr
!k_vecs2d(ik,2)=k_vecs2d(ik,2)*real(nky)/real(nkx)
!k_vecs2d(ik,1)=k_vecs2d(ik,1)*nxx
write(*, "(i5,1x,3f20.10,f7.2)"), ik, k_vecs2d(ik,1),k_vecs2d(ik,2),k_vecs2d(ik,3),k_vecs2d(ik,4)
enddo


 close(unit=100)
 close(unit=101)




!stop

end subroutine set_kvecs_kr_pw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kvecs_kr_pw_nscf (k_vecs2d_irr, nk2d_irr)
real(kind=idp),allocatable, intent(out)	::k_vecs2d_irr(:,:)
integer, intent(out)	::nk2d_irr
integer	::ik,ikn
real	::dummy
real(kind=idp)	:: y_corr, z_corr


if (lysameasx) then
y_corr=1.0d0
else
y_corr=0.9999d0
endif
if (lzsameasx) then
z_corr=1.0d0
else
z_corr=0.9998d0
endif

OPEN(unit=100,file=trim(label)//'/input_k2d_nscf',status='unknown')


write(100,"(i1)") 8
write(100,"(a)") "k2d"
write(100,"(f15.8)") real(nkz)/real(nkx)*z_corr	!to impose that z direction is always different from x-y
write(100,"(f15.8)") real(nky)/real(nkx)*y_corr
write(100,"(3i6)") nkx_plot,nky_plot,1
write(100,"(3i6)") 0,0,0
write(100,"(a1)") "f"

 CLOSE(unit=100)



 call system("kpoints.x <" // trim(label) // "/input_k2d_nscf")
 call system("mv info " // trim(label) // "/info_nscf_2d")
 call system("mv k2d " // trim(label) // "/k2d_nscf" )


OPEN(unit=101,file=trim(label)//'/k2d_nscf',status='unknown')
read(101,"(i)"), nk2d_irr

allocate(k_vecs2d_irr(0:nk2d_irr-1,4))

do ik=0, nk2d_irr-1
!read(101,"(i5,1x,3f20.10,f7.2)"), ikn, k_vecs2d_irr(ikn-1,1),k_vecs2d_irr(ikn-1,2),k_vecs2d_irr(ikn-1,3),k_vecs2d_irr(ikn-1,4)
read(101,*), ikn, k_vecs2d_irr(ikn-1,1),k_vecs2d_irr(ikn-1,2),k_vecs2d_irr(ikn-1,3),k_vecs2d_irr(ikn-1,4)
enddo
 CLOSE(unit=101)

write(*,"(a26,i5)") "Nk points irrid nscf:", nk2d_irr

do ik=0, nk2d_irr-1
k_vecs2d_irr(ik,2)=k_vecs2d_irr(ik,2)*real(nky)/real(nkx)*y_corr
k_vecs2d_irr(ik,3)=k_vecs2d_irr(ik,3)*real(nkz)/real(nkx)*z_corr
k_vecs2d_irr(ik,2)=k_vecs2d_irr(ik,2)*real(nky)/real(nkx)
write(*, "(i5,1x,3f20.10,f7.2)"), ik, k_vecs2d_irr(ik,1),k_vecs2d_irr(ik,2),k_vecs2d_irr(ik,3),k_vecs2d_irr(ik,4)
enddo


 close(unit=100)
 close(unit=101)




!stop

end subroutine set_kvecs_kr_pw_nscf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kr_vecs_nscf(kr_vecs_nscf, nkx_plot, nky_plot)			
real(kind=idp),allocatable,intent(out)	::kr_vecs_nscf(:,:)!,delta_kr_nscf(:)
!complex(kind=idpc),allocatable,intent(out)	::phi_kr_sigma_alpha_nscf(:,:,:,:,:),phi_kr_alpha_sigma_nscf(:,:,:,:,:),kin_energy_kr_nscf(:,:,:,:,:,:),kin_energy2_kr_nscf(:,:,:,:,:,:,:),&
!						kin_energy3_kr_nscf(:,:,:,:,:,:,:),phi_kr_sigma_alpha2_nscf(:,:,:,:,:,:),phi_kr_alpha_sigma2_nscf(:,:,:,:,:,:),&
!						phi_kr_sigma_alpha3_nscf(:,:,:,:,:,:),phi_kr_alpha_sigma3_nscf(:,:,:,:,:,:)
integer::ikx_s,iky_s
integer, intent(in):: nkx_plot, nky_plot

lrednscf=.true.
if (lgreen .or. lborn .or. lg0) lrednscf=.false.

if (.not. lpathnscf .and. lrednscf) then		!reduced set with kpoints.x
 call set_kvecs_kr_pw_nscf(kr_vecs_nscf, nkr_plot) 
! call set_kvecs_kr_pw     (kr_vecs_nscf, nkr_plot, nkx_ploz, nky_plot,2)	!reduced set

 else




 if (.not. lpathnscf) nkr_plot=(nkx_plot)*(nky_plot)	!nscf
! if (lpathnscf) nkr_plot=4*nkpath_plot+2*nkpath_plot_xy+1	!path
 if (lpathnscf) then

 if (lfullbz) then		!ignore reconstruction and build in the original large BZ
nkxpath_plot=nkpath_plot
nkypath_plot=nkpath_plot
nkpath_plot_xy=int(sqrt(2.0D0)*nkpath_plot)
else				!in the reduced BZ
nkxpath_plot=nkpath_plot/nx_kr
nkypath_plot=nkpath_plot/ny_kr
nkpath_plot_xy=int(sqrt(dble(nkxpath_plot**2+nkypath_plot**2)))
endif


if (nx_kr==1 .and. ny_kr==1) nkr_plot=nkxpath_plot+nkypath_plot+nkpath_plot_xy+1	!path GXMG
if (nx_kr==2 .and. ny_kr==1) nkr_plot=2*nkxpath_plot+2*nkypath_plot+1	!path GxmX'g
!if (nx_kr==2 .and. ny_kr==1 .and.       lfullbz) nkr_plot=4*nkxpath_plot+2*nkypath_plot+1	!path GxmX'g
endif

!110 surface
if(surf_type=="110" .and. lpathnscf) then
nkxpath_plot=nkpath_plot
nkypath_plot=int(dble(nkxpath_plot/sqrt(2.0)))
nkr_plot=2*nkxpath_plot+2*nkypath_plot+1
endif
if(surf_type=="210" .and. lpathnscf) then
nkxpath_plot=nkpath_plot
nkypath_plot=int(dble(nkxpath_plot/sqrt(5.0)))
nkr_plot=2*nkxpath_plot+2*nkypath_plot+1
endif
if(surf_type=="111" .and. lpathnscf) then
nkxypath_plot=nkpath_plot
nkypath_plot=int(dble(nkpath_plot*sqrt(3.0)/2.0d0))
nkxpath_plot=nkpath_plot/2
nkr_plot=nkxpath_plot+nkypath_plot+nkxypath_plot+1
endif


print *, "nkr_plot=", nkr_plot

 allocate(kr_vecs_nscf(0:nkr_plot-1,6))
 kr_vecs_nscf=0


 if (.not. lpathnscf) then				!homogeneous mesh
  do ikx=0,nkx_plot-1
   do iky=0,nky_plot-1

   	ik=ikx+(nkx_plot)*iky
!   	ik=iky+(nkx_plot)*ikx
	
if (mod(nkx_plot,2)==0)   kr_vecs_nscf (ik,1)=-0.5d0+real(ikx)/real(nkx_plot)
if (mod(nky_plot,2)==0)   kr_vecs_nscf (ik,2)=-0.5d0+real(iky)/real(nky_plot)

!from -0.5 to 0.5
if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,1)=-0.5d0+real(ikx+0.5d0)/real(nkx_plot)
if (mod(nky_plot,2)==1)   kr_vecs_nscf (ik,2)=-0.5d0+real(iky+0.5d0)/real(nky_plot)

!from 0 to 1
!if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,1)=real(ikx)/real(nkx_plot)	
!if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,2)=real(iky)/real(nkx_plot)

!if (mod(nkx_plot,2)==0)   kr_vecs_nscf (ik,1)=-0.5+real(ikx+0.5)/real(nkx_plot-1)	!k_x 	odd		da cambiare se si vuole la dos
!if (mod(nky_plot,2)==0)   kr_vecs_nscf (ik,2)=-0.5+real(iky+0.5)/real(nky_plot-1)	!k_y

!if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,1)=-0.5+real(ikx)/real(nkx_plot-1)	!k_x	even
!if (mod(nky_plot,2)==1)   kr_vecs_nscf (ik,2)=-0.5+real(iky)/real(nky_plot-1)	!k_y

   enddo
  enddo



!k points in the enlarged brillouin zone -- for 2x1 only (check for even odd for general case)
if (nx_kr>1 .or. ny_kr>1) then
nkr_plot_klarge=nkx_plot*nky_plot*nx_kr*ny_kr
allocate(kr_vecs_nscf_klarge(0:nkr_plot_klarge-1,4))
kr_vecs_nscf_klarge=0
allocate(list_kk(0:nkr_plot_klarge-1,3))
list_kk=0


!!!!!!!nx odd


do ikx=0, nkx_plot*nx_kr-1
do iky=0, nky_plot*ny_kr-1

ik=ikx+nkx_plot*nx_kr*iky

if (mod(nkx_plot*nx_kr,2)==0)   kr_vecs_nscf_klarge(ik,1) =-0.5d0*nx_kr+real(ikx)/real(nkx_plot)
if (mod(nky_plot*ny_kr,2)==0)   kr_vecs_nscf_klarge(ik,2) =-0.5d0*ny_kr+real(iky)/real(nky_plot)

if (mod(nkx_plot*nx_kr,2)==1)   kr_vecs_nscf_klarge(ik,1)=-0.5d0*nx_kr+real(ikx+0.5d0)/real(nkx_plot)
if (mod(nky_plot*ny_kr,2)==1)   kr_vecs_nscf_klarge(ik,2)=-0.5d0*ny_kr+real(iky+0.5d0)/real(nky_plot)

ikx_s=mod(ikx+nkx_plot/2,nkx_plot)	!kx in the reduced BZ
iky_s=mod(iky,nky_plot)			!!kx in the reduced BZ


 list_kk(ik,3)=ikx_s+nkx_plot*iky_s	!index of external momentum	(small)
! list_kk(ik,1)=mod(ikx+nkx_plot*nx_kr-ikx_s-nkx_plot/2,nkx_plot)		!internal momentum x	(large)
! list_kk(ik,2)=mod(iky+nky_plot*ny_kr-iky_s,nky_plot)		!internal momentum y	(large)


enddo
enddo

!do ik=0, nkr_plot_klarge-1
!kr_vecs_nscf_klarge(ik,1)=(kr_vecs_nscf(list_kk(ik,3),1)+list_kk(ik,1))/nx_kr
!kr_vecs_nscf_klarge(ik,2)=(kr_vecs_nscf(list_kk(ik,3),2)+list_kk(ik,2))/ny_kr
!enddo

do ik=0, nkr_plot_klarge-1

!write(*, "(3i5, 2f10.5, i5, 4f10.5)")   ik, list_kk(ik,1),list_kk(ik,2),real(list_kk(ik,1))/nx_kr,real(list_kk(ik,2))/ny_kr,&
!				list_kk(ik,3),kr_vecs_nscf(list_kk(ik,3),1), kr_vecs_nscf(list_kk(ik,3),2),&
!				 kr_vecs_nscf_klarge(ik,1), kr_vecs_nscf_klarge(ik,2) 

enddo


kr_vecs_nscf_klarge(:,4)=1		!weight
kr_vecs_nscf_klarge(:,3)=0		!z



endif

 else

  call set_kr_path(kr_vecs_nscf, nkpath_plot, nkr_plot)

 endif

kr_vecs_nscf(:,4)=1		!weight

 endif

!only for plotting purposes
if (.not. lpathnscf) then 
allocate(kr_vecs_nscf_qpi(0:nkr_plot*nx_kr*ny_kr-1,6))
kr_vecs_nscf_qpi=0

do ikx=0, nkx_plot*nx_kr-1
do iky=0, nky_plot*ny_kr-1

ik=ikx+nkx_plot*nx_kr*iky

if (mod(nkx_plot*nx_kr,2)==0)   kr_vecs_nscf_qpi(ik,1) =-0.5d0+real(ikx)/real(nkx_plot*nx_kr)	!for some reason I need ikx+1
if (mod(nky_plot*ny_kr,2)==0)   kr_vecs_nscf_qpi(ik,2) =-0.5d0+real(iky)/real(nky_plot*ny_kr)

if (mod(nkx_plot*nx_kr,2)==1)   kr_vecs_nscf_qpi(ik,1)=-0.5d0+real(ikx+0.5d0)/real(nkx_plot*nx_kr)
if (mod(nky_plot*ny_kr,2)==1)   kr_vecs_nscf_qpi(ik,2)=-0.5d0+real(iky+0.5d0)/real(nky_plot*ny_kr)


enddo
enddo

kr_vecs_nscf_qpi(:,5)=k_surf(1,1)*kr_vecs_nscf_qpi(:,1)+k_surf(2,1)*kr_vecs_nscf_qpi(:,2)
kr_vecs_nscf_qpi(:,6)=k_surf(1,2)*kr_vecs_nscf_qpi(:,1)+k_surf(2,2)*kr_vecs_nscf_qpi(:,2)

!to reduce to 1st BZ
if(surf_type=="111") then
do ik=0,nkr_plot-1
if (kr_vecs_nscf_qpi(ik,6)+sqrt(3.)*kr_vecs_nscf_qpi(ik,5)-sqrt(2./3.)>=0 .and. (kr_vecs_nscf_qpi(ik,6)>=0)) then
kr_vecs_nscf_qpi(ik,5)=kr_vecs_nscf_qpi(ik,5)-k_surf(1,1)
kr_vecs_nscf_qpi(ik,6)=kr_vecs_nscf_qpi(ik,6)-k_surf(1,2)
endif
if (kr_vecs_nscf_qpi(ik,6)-sqrt(3.)*kr_vecs_nscf_qpi(ik,5)+sqrt(2./3.)<=0 .and. (kr_vecs_nscf_qpi(ik,6)<0)) then
kr_vecs_nscf_qpi(ik,5)=kr_vecs_nscf_qpi(ik,5)+k_surf(2,1)
kr_vecs_nscf_qpi(ik,6)=kr_vecs_nscf_qpi(ik,6)+k_surf(2,2)
endif
if (kr_vecs_nscf_qpi(ik,6)-sqrt(3.)*kr_vecs_nscf_qpi(ik,5)-sqrt(2./3.)>=0 .and. (kr_vecs_nscf_qpi(ik,6)>=0)) then
kr_vecs_nscf_qpi(ik,5)=kr_vecs_nscf_qpi(ik,5)-k_surf(2,1)
kr_vecs_nscf_qpi(ik,6)=kr_vecs_nscf_qpi(ik,6)-k_surf(2,2)
endif
if (kr_vecs_nscf_qpi(ik,6)+sqrt(3.)*kr_vecs_nscf_qpi(ik,5)+sqrt(2./3.)<=0 .and. (kr_vecs_nscf_qpi(ik,6)<0)) then
kr_vecs_nscf_qpi(ik,5)=kr_vecs_nscf_qpi(ik,5)+k_surf(1,1)
kr_vecs_nscf_qpi(ik,6)=kr_vecs_nscf_qpi(ik,6)+k_surf(1,2)
endif
enddo
endif

endif

!if  (.not. lreadham) allocate(phi_kr_sigma_alpha_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1))
!if  (.not. lreadham) allocate(phi_kr_alpha_sigma_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nkr_plot-1))
!if  (.not. lreadham) allocate(phi_kr_sigma_alpha2_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1,n_nn2))
!if  (.not. lreadham) allocate(phi_kr_alpha_sigma2_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nkr_plot-1,n_nn2))
!if  (.not. lreadham) allocate(phi_kr_sigma_alpha3_nscf(0:l_orbital-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1,n_nn3))
!if  (.not. lreadham) allocate(phi_kr_alpha_sigma3_nscf(0:l_orbital_f-1,0:l_orbital-1,0:1,0:1,0:nkr_plot-1,n_nn3))
!if  (.not. lreadham) allocate(kin_energy_kr_nscf		(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1,0:1))
!if  (.not. lreadham) allocate(kin_energy2_kr_nscf		(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1,n_nn2,0:1))
!if  (.not. lreadham) allocate(kin_energy3_kr_nscf		(0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:1,0:nkr_plot-1,n_nn3,0:1))
!if  (.not. lreadham) allocate(delta_kr_nscf(0:nkr_plot-1))

!if  (.not. lreadham) call set_phi_k(phi_kr_sigma_alpha_nscf,phi_kr_alpha_sigma_nscf,kin_energy_kr_nscf, delta_kr_nscf,nkr_plot,kr_vecs_nscf,2)
!if  (.not. lreadham) call set_kin_energy2_kr(kin_energy2_kr_nscf,phi_kr_sigma_alpha2_nscf,phi_kr_alpha_sigma2_nscf,nkr_plot, kr_vecs_nscf)
!if  (.not. lreadham) call set_kin_energy3_kr(kin_energy3_kr_nscf,phi_kr_sigma_alpha3_nscf,phi_kr_alpha_sigma3_nscf,nkr_plot, kr_vecs_nscf)


!correction for the larger full bz
if (lpathnscf .and. lfullbz) kr_vecs_nscf(:,1)=kr_vecs_nscf(:,1)*nx_kr
if (lpathnscf .and. lfullbz) kr_vecs_nscf(:,2)=kr_vecs_nscf(:,2)*ny_kr


kr_vecs_nscf(:,5)=k_surf(1,1)*kr_vecs_nscf(:,1)+k_surf(2,1)*kr_vecs_nscf(:,2)
kr_vecs_nscf(:,6)=k_surf(1,2)*kr_vecs_nscf(:,1)+k_surf(2,2)*kr_vecs_nscf(:,2)

!111 surface
if(surf_type=="111") then
do ik=0,nkr_plot-1
if (kr_vecs_nscf(ik,6)+sqrt(3.)*kr_vecs_nscf(ik,5)-sqrt(2./3.)>=0 .and. (kr_vecs_nscf(ik,6)>=0)) then
kr_vecs_nscf(ik,5)=kr_vecs_nscf(ik,5)-k_surf(1,1)
kr_vecs_nscf(ik,6)=kr_vecs_nscf(ik,6)-k_surf(1,2)
endif
!MK
if (kr_vecs_nscf(ik,6)-sqrt(3.)*kr_vecs_nscf(ik,5)+sqrt(2./3.)<=0 .and. (kr_vecs_nscf(ik,6)<0)) then
kr_vecs_nscf(ik,5)=kr_vecs_nscf(ik,5)+k_surf(2,1)
kr_vecs_nscf(ik,6)=kr_vecs_nscf(ik,6)+k_surf(2,2)
endif
if (kr_vecs_nscf(ik,6)-sqrt(3.)*kr_vecs_nscf(ik,5)-sqrt(2./3.)>=0 .and. (kr_vecs_nscf(ik,6)>=0)) then
kr_vecs_nscf(ik,5)=kr_vecs_nscf(ik,5)-k_surf(2,1)
kr_vecs_nscf(ik,6)=kr_vecs_nscf(ik,6)-k_surf(2,2)
endif
if (kr_vecs_nscf(ik,6)+sqrt(3.)*kr_vecs_nscf(ik,5)+sqrt(2./3.)<=0 .and. (kr_vecs_nscf(ik,6)<0)) then
kr_vecs_nscf(ik,5)=kr_vecs_nscf(ik,5)+k_surf(1,1)
kr_vecs_nscf(ik,6)=kr_vecs_nscf(ik,6)+k_surf(1,2)
endif

enddo
endif


end subroutine set_kr_vecs_nscf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine do_kr_nscf(b_kr,lambda_kr, mu_kr)
real(kind=idp)	:: b_kr(0:n_sites_kr-1),lambda_kr(0:n_sites_kr-1), mu_kr
integer:: ind1,indp, ind, indp2, jmax,imax
real(kind=idp):: emax, emin
 CHARACTER(5) :: ak

 do ix=0,nx_kr-1	!x
   do iy=0,ny_kr-1	!y
    do iz=0,nz-1	!z
 ind=index_xyz(ix,iy,iz,nx_kr,ny_kr)
!if (iz==0 .or. iz==nz-1) b_kr(ind)=b_kr0
!if (iz==1 .or. iz==nz-2) b_kr(ind)=b_kr1
   enddo
  enddo
 enddo


  write(*,"(a30)") ""
  write(*,"(a30)") "Doing NSCF with:"
  write(*,"(a10,f15.8)") "b, lambda ="
do ix=0,nx_kr-1
 do iy=0,ny_kr-1
  do iz=0,nz-1
   i=index_xyz(ix,iy,iz,nx_kr,ny_kr)
  write(*,"(3i5,2f15.8)") ix,iy,iz,b_kr(i), lambda_kr(i)
  enddo
enddo
enddo
  write(*,"(a10,f15.8)") "mu=",mu_kr
  write(*,"(a10,i15)") "k points=",nkr_plot
  
  
!size_h_kr_red_e=min(size_h_kr_red_e,size_h_kr)

size_h_kr_red_e=size_h_kr


if (lpathnscf)  size_h_kr_red_e=size_h_kr
!if (size_h_kr_red_e==size_h_kr) emin=-1e5			!no reduction in energy: I keep all states
!if (size_h_kr_red_e.ne.size_h_kr) emin=omega_min-10*T_g	!reduction in energy: I start from the minimum energy needed

if (size_h_kr_red_e.ne.size_h_kr) print *, "WARNING: energy reduction"
if (size_h_kr_red_e==size_h_kr) print *, "No energy reduction"


write(*, "(a50,4i10)"), "Size W_kr_tot_nscf (basis,eigenvalues,k pts, total):",size_h_kr_red,size_h_kr_red_e,nkr_plot, size_h_kr_red*size_h_kr_red_e*nkr_plot
!write(*, "(a10,f15.5)"), "emin=", emin

!  allocate(corr_fct_kr_nscf(0:size_h_kr-1,0:size_h_kr-1))
		allocate(e_kr_tot_nscf(0:size_h_kr_red_e-1,0:nkr_plot-1))
		allocate(w_kr_tot_nscf(0:size_h_kr_red-1,0:size_h_kr_red_e-1,0:nkr_plot-1))
if (lpathnscf)  allocate(ew_kr_tot_nscf(0:size_h_kr_red_e-1,0:nkr_plot-1,1:3,0:(l_orb_tot)-1))
if (lpathnscf)  allocate(ew_kr_tot_nscf_rec(0:size_h_kr_red_e-1,0:nkr_plot-1,1:5,0:(l_orb_tot)-1))
if (lpathnscf)  allocate(spin_kr_tot_nscf(0:size_h_kr_red_e-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)  allocate(spin_kr_tot_nscf_rot(0:size_h_kr_red_e-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)  allocate(spin_kr_tot_p_nscf(0:size_h_kr-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)  allocate(spin_kr_tot_m_nscf(0:size_h_kr-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)  allocate(spin_kr_tot_p_nscf_rot(0:size_h_kr-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)  allocate(spin_kr_tot_m_nscf_rot(0:size_h_kr-1,0:nkr_plot-1,1:3,0:1,0:nz_plot-1))
if (lpathnscf)	allocate(w_kr_tot_nscf_pm(0:size_h_kr/2-1,0:size_h_kr/2-1,0:nkr_plot-1,2))
if (lpathnscf)	allocate(e_kr_tot_nscf_pm(0:size_h_kr/2-1,0:nkr_plot-1,2))
if (lpathnscf)  allocate(ew_kr_tot_p_nscf(0:size_h_kr-1,0:nkr_plot-1,1:3,0:(l_orb_tot)-1))
if (lpathnscf)  allocate(ew_kr_tot_m_nscf(0:size_h_kr-1,0:nkr_plot-1,1:3,0:(l_orb_tot)-1))
  !allocate(H_kr_tot_nscf(0:size_h_kr-1,0:size_h_kr-1,0:nkr_plot-1))
!  allocate(v_kr_nscf(0:size_h_kr-1,0:size_h_kr-1))
!  allocate(e_corr_fct_kr_nscf(0:size_h_kr-1,0:nkr_plot-1))
  
 
  
!  occup_kr=0
!  fermien_kr(:,1)=100
!  fermien_kr(:,2)=-100
 ! e_corr_fct_kr_nscf=0

e_kr_tot_nscf=0
w_kr_tot_nscf=0
if (lpathnscf) e_kr_tot_nscf_pm=0
if (lpathnscf) w_kr_tot_nscf_pm=0
if (lpathnscf) ew_kr_tot_nscf=0
if (lpathnscf) ew_kr_tot_nscf_rec=0
jmax=0
emax=1e4



!Symmetries
allocate(M_yz(0:dim_hk-1,0:dim_hk-1))
allocate(w_yz(0:dim_hk-1,0:dim_hk-1))
allocate(e_yz(0:dim_hk-1))
allocate(w_yz_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_yz_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_yz_m_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_yz_p_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(M_xy(0:dim_hk-1,0:dim_hk-1))
allocate(w_xy(0:dim_hk-1,0:dim_hk-1))
allocate(e_xy(0:dim_hk-1))
allocate(M_x(0:dim_hk-1,0:dim_hk-1))
allocate(w_x(0:dim_hk-1,0:dim_hk-1))
allocate(e_x(0:dim_hk-1))
allocate(M_y(0:dim_hk-1,0:dim_hk-1))
allocate(w_y(0:dim_hk-1,0:dim_hk-1))
allocate(e_y(0:dim_hk-1))
allocate(w_xy_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_xy_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_y_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_y_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_x_m(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_x_p(0:dim_hk-1,0:dim_hk/2-1))
allocate(w_xy_m_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_xy_p_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_y_m_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_y_p_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_x_m_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_x_p_kr(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(ham_kr_tmp_p(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(ham_kr_tmp_m(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(ham_kr_tmp_pm(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(ham_kr_tmp_mp(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(ham_kr_p(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(ham_kr_m(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(ham_kr_pm(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(ham_kr_mp(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(w_kr_p(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_kr_m(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_kr_p_tmp(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(w_kr_m_tmp(0:size_h_kr/2-1,0:size_h_kr/2-1))
allocate(e_kr_p(0:size_h_kr/2-1))
allocate(e_kr_m(0:size_h_kr/2-1))
allocate(w_x_p_kr_pi(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(w_x_m_kr_pi(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(e_x_kr(0:size_h_kr-1))
allocate(e_x_kr_pi(0:size_h_kr-1))
allocate(w_x_kr(0:size_h_kr-1,0:size_h_kr-1))
allocate(w_x_kr_pi(0:size_h_kr-1,0:size_h_kr-1))
allocate(u_symm_kr_p(0:size_h_kr-1,0:size_h_kr/2-1))
allocate(u_symm_kr_m(0:size_h_kr-1,0:size_h_kr/2-1))



M_x(:,:)=simm_oper(:,:,3)		!3:M_x	4:M_y 5:M_z
M_y(:,:)=simm_oper(:,:,4)		!3:M_x	4:M_y 5:M_z
M_xy(:,:)=simm_oper(:,:,6)		!6:M_xy
M_yz(:,:)=simm_oper(:,:,5)		!5:M_z
w_x=0
w_y=0
w_xy=0
w_yz=0
call DIAGONALIZE( 'V', dim_hk, dim_hk, M_x, e_x, w_x)!reflections: hermitian operator (times i, that I have removed)
call DIAGONALIZE( 'V', dim_hk, dim_hk, M_y, e_y, w_y)!reflections: hermitian operator (times i, that I have removed)
call DIAGONALIZE( 'V', dim_hk, dim_hk, M_xy, e_xy, w_xy)!reflections: hermitian operator (times i, that I have removed)
call DIAGONALIZE( 'V', dim_hk, dim_hk, M_yz, e_yz, w_yz)!reflections: hermitian operator (times i, that I have removed)
do i=0, dim_hk/2-1
w_xy_m(:,i)=w_xy(:,i)		!eigenstates -
w_xy_p(:,i)=w_xy(:,i+dim_hk/2)	!eigenstates +
w_x_m(:,i)=w_x(:,i)		!eigenstates -
w_x_p(:,i)=w_x(:,i+dim_hk/2)	!eigenstates +
w_y_m(:,i)=w_y(:,i)		!eigenstates -
w_y_p(:,i)=w_y(:,i+dim_hk/2)	!eigenstates +
w_yz_m(:,i)=w_yz(:,i)		!eigenstates -
w_yz_p(:,i)=w_yz(:,i+dim_hk/2)	!eigenstates +
enddo

!w_x_p_kr=0
w_y_p_kr=0
w_xy_p_kr=0
!w_x_m_kr=0
w_y_m_kr=0
w_xy_m_kr=0
w_yz_m_kr=0
w_yz_m_kr=0
do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   ind1= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz)   	
   ind2= index_cso(icf,is,o1)
   do j=0, dim_hk/2-1
   ind3= index_xyzcso2(ix,iy,iz,j,nx_kr,ny_kr,nz)
!   ind4= index_xyzcso2(1-ix,iy,iz,j,nx_kr,ny_kr,nz)
!w_x_p_kr(ind1,ind3)=w_x_p(ind2,j)*exp(2*pi*(0,1)*ix*kr_vecs_nscf(ik,1))		!comes from reflection operator
!w_x_m_kr(ind1,ind3)=w_x_m(ind2,j)*exp(2*pi*(0,1)*ix*kr_vecs_nscf(ik,1))
w_y_p_kr(ind1,ind3)=w_y_p(ind2,j)
w_y_m_kr(ind1,ind3)=w_y_m(ind2,j)
w_xy_p_kr(ind1,ind3)=w_xy_p(ind2,j)
w_xy_m_kr(ind1,ind3)=w_xy_m(ind2,j)
w_yz_p_kr(ind1,ind3)=w_yz_p(ind2,j)
w_yz_m_kr(ind1,ind3)=w_yz_m(ind2,j)
  enddo
enddo
enddo
enddo
enddo
enddo
enddo


allocate(M_x_kr(0:size_h_kr-1,0:size_h_kr-1))
allocate(M_x_kr_pi(0:size_h_kr-1,0:size_h_kr-1))
M_x_kr=0
M_x_kr_pi=0

do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   ind1= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz)   	
   ind2= index_cso(icf,is,o1)
ix2=ix
iy2=iy
iz2=iz
 do icf2=0,1
  do is2=0,1
   do o2=0,l_orb(icf)-1
   ind3= index_xyzcso(ix2,iy2,iz2,icf2,is2,o2,nx_kr,ny_kr,nz)   	
   ind4= index_cso(icf2,is2,o2)

M_x_kr(ind1,ind3)=M_x(ind2,ind4)
M_x_kr(ind1,ind3)=m_x(ind2,ind4)
if (ix==0) M_x_kr_pi(ind1,ind3)=+m_x(ind2,ind4)
if (ix==1) M_x_kr_pi(ind1,ind3)=-m_x(ind2,ind4)


enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


call DIAGONALIZE( 'V', size_h_kr, size_h_kr, M_x_kr, e_x_kr, w_x_kr)		!reflections: hermitian operator (times i, that I have removed)
call DIAGONALIZE( 'V', size_h_kr, size_h_kr, M_x_kr_pi, e_x_kr_pi, w_x_kr_pi)	!reflections: hermitian operator (times i, that I have removed)


do i=0, size_h_kr/2-1
w_x_m_kr(:,i)=w_x_kr(:,i)		!eigenstates -
w_x_p_kr(:,i)=w_x_kr(:,i+size_h_kr/2)		!eigenstates -
w_x_m_kr_pi(:,i)=w_x_kr_pi(:,i)		!eigenstates -
w_x_p_kr_pi(:,i)=w_x_kr_pi(:,i+size_h_kr/2)		!eigenstates -
enddo




!print *, "u-1_kr * u_kr"
! call write_mat(ham_kr_p,size_h_kr/2)
! call write_mat(ham_kr_m,size_h_kr/2)



do ik=0, nkr_plot-1
if (mod(ik,1000)==0) write(*, "(a5,i10,f15.5,i5)") "k=", ik, emax, jmax

if (lpathnscf) then
ham_kr_p=0
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, w_x_p_kr, w_x_p_kr, ham_kr_p)
ham_kr_m=0
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, w_x_m_kr, w_x_m_kr, ham_kr_m)
ham_kr_pm=0
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, w_x_p_kr, w_x_m_kr, ham_kr_pm)
ham_kr_mp=0
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, w_x_m_kr, w_x_p_kr, ham_kr_mp)

do i=0, size_h_kr/2-1
if (abs(ham_kr_p(i,i)-1)>1e-5) print *, "u-1u +", i, ham_kr_p(i,i)
if (abs(ham_kr_m(i,i)-1)>1e-5) print *, "u-1u -", i, ham_kr_m(i,i)
do j=0, size_h_kr/2-1
if (j.ne.i .and. abs(ham_kr_p(i,j))>1e-5) print *, "u-1u ++", i, ham_kr_p(i,j)
if (j.ne.i .and. abs(ham_kr_m(i,j))>1e-5) print *, "u-1u --", i, ham_kr_m(i,j)
if (abs(ham_kr_pm(i,j))>1e-5) print *, "u-1u +-", i, ham_kr_p(i,j)
if (abs(ham_kr_mp(i,j))>1e-5) print *, "u-1u -+", i, ham_kr_m(i,j)
enddo
enddo
endif

!   write (*, "(i5, 4f15.8)"), ik, kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2),kr_vecs_nscf(ik,3),kr_vecs_nscf(ik,4)
!print *, ik
!!!build and diag hamiltionian

!if (.not. lreadham .and. .not. l2dmodel)  call build_ham_kr(ik,Ham_kr,b_kr,lambda_kr, mu_kr,.true.,nkr_plot,phi_kr_sigma_alpha_nscf,phi_kr_alpha_sigma_nscf, kin_energy_kr_nscf,&
! kin_energy2_kr_nscf,kin_energy3_kr_nscf,phi_kr_sigma_alpha2_nscf,phi_kr_alpha_sigma2_nscf,phi_kr_sigma_alpha3_nscf,phi_kr_alpha_sigma3_nscf,kr_vecs_nscf)
!if (lreadham .and. .not. l2dmodel)        call build_ham_kr_read(ik,Ham_kr,b_kr,lambda_kr, mu_kr,.true.,nkr_plot,kr_vecs_nscf)
if (l2dmodel)				 call build_ham_kr_2d(ik,Ham_kr,nkr_plot,kr_vecs_nscf)
if (.not. l2dmodel) 			 call build_ham_kr2(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.)

if (lpathnscf) then
ham_kr_tmp_p=0
ham_kr_tmp_m=0
ham_kr_tmp_pm=0
ham_kr_tmp_mp=0
ham_kr_p=0
ham_kr_m=0
ham_kr_pm=0
ham_kr_mp=0
u_symm_kr_p=0
u_symm_kr_m=0

!choose symmetry
if (abs(kr_vecs_nscf(ik,1))<1e-5  .and. surf_type=="001") then			!x=0
u_symm_kr_p=w_x_p_kr
u_symm_kr_m=w_x_m_kr
elseif (abs((kr_vecs_nscf(ik,1)-0.5d0))<1e-5) then			!x=0.5
u_symm_kr_p=w_x_p_kr_pi
u_symm_kr_m=w_x_m_kr_pi
elseif ((abs(kr_vecs_nscf(ik,2))<1e-5 .or. abs(kr_vecs_nscf(ik,2)-0.5)<1e-5) .and. surf_type=="001") then		!y=0 or y=0.5
u_symm_kr_p=w_y_p_kr
u_symm_kr_m=w_y_m_kr
elseif (abs(kr_vecs_nscf(ik,1)-kr_vecs_nscf(ik,2))<1e-5 .and. nx_kr==1 .and. ny_kr==1  .and. surf_type=="001") then			!y=x
u_symm_kr_p=w_xy_p_kr
u_symm_kr_m=w_xy_m_kr
endif
!110 210
if ((abs(kr_vecs_nscf(ik,1))<1e-5 .or. abs(kr_vecs_nscf(ik,1)-0.5)<1e-5) .and. (surf_type=="110" .or. surf_type=="210")) then	!bar k_x=k_x
u_symm_kr_p=w_yz_p_kr	!it is actually m_z
u_symm_kr_m=w_yz_m_kr
endif
if (abs(kr_vecs_nscf(ik,2))<1e-5 .and. surf_type=="110") then	!bar k_y=ky-kx
u_symm_kr_p=w_xy_p_kr
u_symm_kr_m=w_xy_m_kr
endif
if (abs(kr_vecs_nscf(ik,2)-0.5)<1e-5 .and. surf_type=="110") then	!bar k_y=ky-kx
u_symm_kr_p=0
u_symm_kr_m=0
endif
!111
if (abs(kr_vecs_nscf(ik,1)-kr_vecs_nscf(ik,2))<1e-5 .and. surf_type=="111") then	!bar k_y=ky-kx
u_symm_kr_p=w_xy_p_kr
u_symm_kr_m=w_xy_m_kr
endif
if (abs(kr_vecs_nscf(ik,1)-kr_vecs_nscf(ik,2))>1e-5 .and. surf_type=="111") then	!bar k_y=ky-kx
u_symm_kr_p=0
u_symm_kr_m=0
endif



call matmul2 (size_h_kr,size_h_kr/2,size_h_kr, ham_kr,u_symm_kr_p, ham_kr_tmp_p) 	
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr,u_symm_kr_p, ham_kr_tmp_p, ham_kr_p)
call matmul2 (size_h_kr,size_h_kr/2,size_h_kr, ham_kr,u_symm_kr_m, ham_kr_tmp_m) 	
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, u_symm_kr_m, ham_kr_tmp_m, ham_kr_m) 	
call matmul2 (size_h_kr,size_h_kr/2,size_h_kr, ham_kr,u_symm_kr_m, ham_kr_tmp_pm) 	
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, u_symm_kr_p, ham_kr_tmp_pm, ham_kr_pm) 	
call matmul2 (size_h_kr,size_h_kr/2,size_h_kr, ham_kr,u_symm_kr_p, ham_kr_tmp_mp) 	
call matmul2hc2 (size_h_kr/2,size_h_kr/2,size_h_kr, u_symm_kr_m, ham_kr_tmp_mp, ham_kr_mp) 	



do i=0, size_h_kr/2-1
do j=0, size_h_kr/2-1
!if (abs(ham_kr_pm(i,j))>1e-5) print *, "warning, H_kr +-",ik, ham_kr_pm(i,j)
!if (abs(ham_kr_mp(i,j))>1e-5) print *, "warning, H_kr -+",ik, ham_kr_mp(i,j)
enddo
enddo
endif

!Diag whole Ham
  e_kr=0
  call diagonalize_ham_kr(Ham_kr,e_kr,w_kr)

!Diag Ham+ Ham - ...
if (lpathnscf) then
  e_kr_p=0
  e_kr_m=0
  call diagonalize_ham_kr2(Ham_kr_p,e_kr_p,w_kr_p_tmp)
  call diagonalize_ham_kr2(Ham_kr_m,e_kr_m,w_kr_m_tmp)

!...and back to the original basis (I need it to know the orbital and z weight)
  w_kr_p=0
  w_kr_m=0
call matmul2 (size_h_kr,size_h_kr/2,size_h_kr/2, u_symm_kr_p, w_kr_p_tmp, w_kr_p) 	
call matmul2 (size_h_kr,size_h_kr/2,size_h_kr/2, u_symm_kr_m, w_kr_m_tmp, w_kr_m) 	
endif


 if (.not. lscfkr) mu_new_kr=mu_kr
 
!PATH/NSCF nored

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (size_h_kr_red_e==size_h_kr) then
  
  e_kr_tot_nscf(:,ik)=e_kr(:)
if (lpathnscf)  e_kr_tot_nscf_pm(:,ik,1)=e_kr_p(:)
if (lpathnscf)  e_kr_tot_nscf_pm(:,ik,2)=e_kr_m(:)

do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz_g0-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   ind1= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz)   	!ind of w_kr -->nz
   indp= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)	!ind of w_kr_tot -->nz_g0

  w_kr_tot_nscf(indp,:,ik)=w_kr(ind1,:)
  
enddo
enddo
enddo
enddo
enddo
enddo


!!!!!!!!!!!!!!!!!
!I fourier transform to go into the large BZ
if (lrec) then
w_kr_rec=0
do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   ind1= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz)   	!ind of w_kr -->nz
   ind=index_co(icf,o1)
   ind2= index_zcso(iz,icf,is,o1,nz)   	!ind of w_kr -->nz

do i=0, size_h_kr-1

w_kr_rec(ind2,i)=w_kr_rec(ind2,i)+w_kr(ind1,i)*exp(2*pi*(0,1)*kr_vecs_nscf(ik,1)*ix/nx_kr)	

!dos_arpes_klarge(ik,icf,iz,o1)=dos_arpes_klarge(ik,icf,iz,o1)-1/pi*aimag(G0_k(indp,indp2,list_kk(ik,3))*&
!			exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))


enddo
enddo
enddo
enddo
enddo
enddo
enddo
endif




!weight on each orbital, and on first layer
if (lpathnscf .and. .not. l2dmodel) then
if (ik==0) ew_kr_tot_nscf=0
if (ik==0) ew_kr_tot_nscf_rec=0
if (ik==0) ew_kr_tot_p_nscf=0
if (ik==0) ew_kr_tot_m_nscf=0
if (ik==0) spin_kr_tot_nscf=0
if (ik==0) spin_kr_tot_nscf_rot=0
if (ik==0) spin_kr_tot_p_nscf=0
if (ik==0) spin_kr_tot_m_nscf=0
if (ik==0) spin_kr_tot_p_nscf_rot=0
if (ik==0) spin_kr_tot_m_nscf_rot=0



do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   ind1= index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz)   	!ind of w_kr -->nz
   ind=index_co(icf,o1)
   ind3= index_zcso(iz,icf,is,o1,nz)   	!ind of w_kr -->nz

do i=0, size_h_kr-1

 				ew_kr_tot_nscf(i,ik,1,ind)=ew_kr_tot_nscf(i,ik,1,ind)+abs(w_kr(ind1,i))**2	!total weight on each orbital
 if (iz==0 .or. iz==nz-1)	ew_kr_tot_nscf(i,ik,2,ind)=ew_kr_tot_nscf(i,ik,2,ind)+abs(w_kr(ind1,i))**2/2	!weight on the first layer on each orbital; first layer weight:  (first+last)/2
 if (iz==1 .or. iz==nz-2)	ew_kr_tot_nscf(i,ik,3,ind)=ew_kr_tot_nscf(i,ik,3,ind)+abs(w_kr(ind1,i))**2/2	!weight on the second layer on each orbital; first+1 layer weight:  (first+1+last-1)/2


!total weight on each orbital in the reconstructed scheme
if (lrec)				ew_kr_tot_nscf_rec(i,ik,1,ind)=ew_kr_tot_nscf_rec(i,ik,1,ind)+abs(w_kr_rec(ind3,i))**2	
if (lrec .and. (iz==0 .or. iz==nz-1))	ew_kr_tot_nscf_rec(i,ik,2,ind)=ew_kr_tot_nscf_rec(i,ik,2,ind)+abs(w_kr_rec(ind3,i))**2/2	!weight on the first layer on each orbital
if (lrec .and. (iz==1 .or. iz==nz-2))	ew_kr_tot_nscf_rec(i,ik,3,ind)=ew_kr_tot_nscf_rec(i,ik,3,ind)+abs(w_kr_rec(ind3,i))**2/2	!weight on the second layer on each orbital
if (lrec .and. (iz==2 .or. iz==nz-3))	ew_kr_tot_nscf_rec(i,ik,4,ind)=ew_kr_tot_nscf_rec(i,ik,4,ind)+abs(w_kr_rec(ind3,i))**2/2	!weight on the second layer on each orbital
if (lrec .and. (iz==3 .or. iz==nz-4))	ew_kr_tot_nscf_rec(i,ik,5,ind)=ew_kr_tot_nscf_rec(i,ik,5,ind)+abs(w_kr_rec(ind3,i))**2/2	!weight on the second layer on each orbital


do ixyz=1,3
  do is2=0,1
   do o2=0,l_orb(icf)-1
   ind2= index_xyzcso(ix,iy,iz,icf,is2,o2,nx_kr,ny_kr,nz)   	!ind of w_kr -->nz

if (iz<nz_plot) spin_kr_tot_nscf(i,ik,ixyz,icf,iz)=spin_kr_tot_nscf(i,ik,ixyz,icf,iz)+conjg(w_kr(ind2,i))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)*w_kr(ind1,i)
!spin_arpes(ik,icf,ixyz,iz,0)=spin_arpes(ik,icf,ixyz,iz,0)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) 

enddo
enddo
enddo	!ixyz

enddo	!i



do i=0, size_h_kr/2-1
 		ew_kr_tot_p_nscf(i,ik,1,ind)=ew_kr_tot_p_nscf(i,ik,1,ind)+abs(w_kr_p(ind1,i))**2	!total weight on each orbital
 if (iz==0)	ew_kr_tot_p_nscf(i,ik,2,ind)=ew_kr_tot_p_nscf(i,ik,2,ind)+abs(w_kr_p(ind1,i))**2/2	!weight on the first layer on each orbital
 if (iz==1)	ew_kr_tot_p_nscf(i,ik,3,ind)=ew_kr_tot_p_nscf(i,ik,3,ind)+abs(w_kr_p(ind1,i))**2/2	!weight on the second layer on each orbital

 		ew_kr_tot_m_nscf(i,ik,1,ind)=ew_kr_tot_m_nscf(i,ik,1,ind)+abs(w_kr_m(ind1,i))**2	!total weight on each orbital
 if (iz==0)	ew_kr_tot_m_nscf(i,ik,2,ind)=ew_kr_tot_m_nscf(i,ik,2,ind)+abs(w_kr_m(ind1,i))**2/2	!weight on the first layer on each orbital
 if (iz==1)	ew_kr_tot_m_nscf(i,ik,3,ind)=ew_kr_tot_m_nscf(i,ik,3,ind)+abs(w_kr_m(ind1,i))**2/2	!weight on the second layer on each orbital



do ixyz=1,3
  do is2=0,1
   do o2=0,l_orb(icf)-1
   ind2= index_xyzcso(ix,iy,iz,icf,is2,o2,nx_kr,ny_kr,nz)   	!ind of w_kr -->nz

if (iz<nz_plot) spin_kr_tot_p_nscf(i,ik,ixyz,icf,iz)=spin_kr_tot_p_nscf(i,ik,ixyz,icf,iz)+conjg(w_kr_p(ind2,i))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)*w_kr_p(ind1,i)
if (iz<nz_plot) spin_kr_tot_m_nscf(i,ik,ixyz,icf,iz)=spin_kr_tot_m_nscf(i,ik,ixyz,icf,iz)+conjg(w_kr_m(ind2,i))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)*w_kr_m(ind1,i)
!spin_arpes(ik,icf,ixyz,iz,0)=spin_arpes(ik,icf,ixyz,iz,0)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) 

enddo
enddo
enddo	!ixyz


enddo	!i
  
enddo
enddo
enddo
enddo
enddo
enddo

!rotates for surface orientation
do i=0, size_h_kr-1
do iz=0, nz_plot-1
do ixyz=1,3
do icf=0,1
do ixyz2=1,3
spin_kr_tot_nscf_rot(i,ik,ixyz,icf,iz)=spin_kr_tot_nscf_rot(i,ik,ixyz,icf,iz)+rot_mat_surf(ixyz,ixyz2)*spin_kr_tot_nscf(i,ik,ixyz2,icf,iz)
spin_kr_tot_p_nscf_rot(i,ik,ixyz,icf,iz)=spin_kr_tot_p_nscf_rot(i,ik,ixyz,icf,iz)+rot_mat_surf(ixyz,ixyz2)*spin_kr_tot_p_nscf(i,ik,ixyz2,icf,iz)
spin_kr_tot_m_nscf_rot(i,ik,ixyz,icf,iz)=spin_kr_tot_m_nscf_rot(i,ik,ixyz,icf,iz)+rot_mat_surf(ixyz,ixyz2)*spin_kr_tot_m_nscf(i,ik,ixyz2,icf,iz)
enddo
enddo
enddo
enddo
enddo

endif





!NSCF with cutting procedure for green
!else 

!j=0 
!imax=size_h_kr-1
!do i=0, size_h_kr-1

!if(e_kr(i)>emin .and. j<size_h_kr_red_e) then 
 
!  e_kr_tot_nscf(j,ik)=e_kr(i)
  
!do iz=0, nz_g0-1
! do icf=0,1
!  do is=0,1
!   do o1=0,l_orb(icf)-1
!   ind1= index_zcso(iz,icf,is,o1,nz)   	!ind of w_kr -->nz
!   indp= index_zcso(iz,icf,is,o1,nz_g0)	!ind of w_kr_tot -->nz_scf

!  w_kr_tot_nscf(indp,j,ik)=w_kr(ind1,i)
  
! if(j==size_h_kr_red_e-1) imax=i+1
 
!enddo
!enddo
!enddo
!enddo

!j=j+1
!endif
!enddo

!if (e_kr_tot_nscf(size_h_kr_red_e-1,ik)<emax) print *, e_kr_tot_nscf(:,ik)
!if (e_kr_tot_nscf(size_h_kr_red_e-1,ik)<emax) print *, e_kr(:)
!if (e_kr_tot_nscf(size_h_kr_red_e-1,ik)<emax) emax=e_kr_tot_nscf(size_h_kr_red_e-1,ik)
!if (e_kr(min(imax,size_h_kr-1))<emax) emax=e_kr(min(imax,size_h_kr-1))			!the first eigenvalue left out (i+1)
!if (j>jmax) jmax=j

!endif !path



enddo	! k points

print *, "jmax=", jmax
OPEN(unit=100,file=trim(label)//'/emax',status='unknown')
write (*, "(a10,f15.5)"), "emax=", emax
write (100, "(a10,f15.5)"), "emax=", emax
 close (unit=100)

!!
!if (lgreen .or. lborn .or. lg0) call w_kr_analysis(e_kr_tot_nscf,w_kr_tot_nscf)

!!!compute spectral weights
!if (lspectrfunct_kr .and. nz==nz_scf) call spectralfunction_kr_nscf(e_kr_tot_nscf,w_kr_tot_nscf)

!write bands and corr
if (lpathnscf) then
if (lpathnscf .and. (.not. lrec .or. .not. lfullbz)) then
OPEN(unit=28,file=trim(label)//'/bands_kr_nscf',status='unknown')

 write(28,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr-1
     do ik= 0, nkr_plot-1
write(28,"(i5,3f10.5, f17.10,15e12.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf(i,ik),&
						(ew_kr_tot_nscf(i,ik,1,j),j=0,l_orb_tot-1),(ew_kr_tot_nscf(i,ik,2,j), j=0,l_orb_tot-1),(ew_kr_tot_nscf(i,ik,3,j), j=0,l_orb_tot-1)
     enddo
write(28,"(a10)") "   "
write(28,"(a10)") "   "
    enddo
 close(unit=28)
endif

print *, lpathnscf, lrec, lfullbz
if(lpathnscf .and. lrec .and. lfullbz) then
OPEN(unit=28,file=trim(label)//'/bands_kr_nscf_rec',status='unknown')

 write(28,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr-1
     do ik= 0, nkr_plot-1
write(28,"(i5,3f10.5, f17.10,40e13.5)") ik, kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), kr_vecs_nscf(ik,3),e_kr_tot_nscf(i,ik),&
						(ew_kr_tot_nscf_rec(i,ik,1,j),j=0,l_orb_tot-1),&
						(ew_kr_tot_nscf_rec(i,ik,2,j),j=0,l_orb_tot-1),&
						(ew_kr_tot_nscf_rec(i,ik,3,j),j=0,l_orb_tot-1),&
						(ew_kr_tot_nscf_rec(i,ik,4,j),j=0,l_orb_tot-1),&
						(ew_kr_tot_nscf_rec(i,ik,5,j),j=0,l_orb_tot-1)
     enddo
write(28,"(a10)") "   "
write(28,"(a10)") "   "
    enddo
 close(unit=28)
endif


OPEN(unit=227,file=trim(label)//'/spins_kr_nscf',status='unknown')
 write(227,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr-1
     do ik= 0, nkr_plot-1
write(227,"(i5,3f10.5, f17.10,30e13.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf(i,ik),&
						sum(spin_kr_tot_nscf_rot(i,ik,1,0,:)), sum(spin_kr_tot_nscf_rot(i,ik,2,0,:)),sum(spin_kr_tot_nscf_rot(i,ik,3,0,:)),&	!cx,cy,cz summed over z
						sum(spin_kr_tot_nscf_rot(i,ik,1,1,:)), sum(spin_kr_tot_nscf_rot(i,ik,2,1,:)),sum(spin_kr_tot_nscf_rot(i,ik,3,1,:)),&	!fx,fy,fz summed over z
						(spin_kr_tot_nscf_rot(i,ik,1,0,0)), (spin_kr_tot_nscf_rot(i,ik,2,0,0)),(spin_kr_tot_nscf_rot(i,ik,3,0,0)),&	!cx,cy,cz z=0
						(spin_kr_tot_nscf_rot(i,ik,1,1,0)), (spin_kr_tot_nscf_rot(i,ik,2,1,0)),(spin_kr_tot_nscf_rot(i,ik,3,1,0))	!fx,fy,fz z=0
						
     enddo
write(227,"(a10)") "   "
write(227,"(a10)") "   "
    enddo
 close(unit=227)
 




OPEN(unit=128,file=trim(label)//'/bands_kr_nscf_p',status='unknown')
OPEN(unit=129,file=trim(label)//'/bands_kr_nscf_m',status='unknown')
write(128,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr/2-1
     do ik= 0, nkr_plot-1
      write(128,"(i5,3f10.5, f17.10,15e12.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf_pm(i,ik,1),& 
						(ew_kr_tot_p_nscf(i,ik,1,j),j=0,l_orb_tot-1),(ew_kr_tot_p_nscf(i,ik,2,j), j=0,l_orb_tot-1),(ew_kr_tot_p_nscf(i,ik,3,j), j=0,l_orb_tot-1)
      write(129,"(i5,3f10.5, f17.10,15e12.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf_pm(i,ik,2),& 
						(ew_kr_tot_m_nscf(i,ik,1,j),j=0,l_orb_tot-1),(ew_kr_tot_m_nscf(i,ik,2,j), j=0,l_orb_tot-1),(ew_kr_tot_m_nscf(i,ik,3,j), j=0,l_orb_tot-1)

     enddo
write(128,"(a10)") "   "
write(128,"(a10)") "   "
write(129,"(a10)") "   "
write(129,"(a10)") "   "
enddo

 close(unit=128)
 close(unit=129)

OPEN(unit=228,file=trim(label)//'/spins_kr_nscf_p',status='unknown')
OPEN(unit=229,file=trim(label)//'/spins_kr_nscf_m',status='unknown')
 write(228,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
 write(229,"(3a10, 5a15)") "#k_x","k_y","k_z", "E"
    do i=0,size_h_kr/2-1
     do ik= 0, nkr_plot-1
write(228,"(i5,3f10.5, f17.10,30e13.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf_pm(i,ik,1),&
						sum(spin_kr_tot_p_nscf_rot(i,ik,1,0,:)), sum(spin_kr_tot_p_nscf_rot(i,ik,2,0,:)),sum(spin_kr_tot_p_nscf_rot(i,ik,3,0,:)),&	!cx,cy,cz summed over z
						sum(spin_kr_tot_p_nscf_rot(i,ik,1,1,:)), sum(spin_kr_tot_p_nscf_rot(i,ik,2,1,:)),sum(spin_kr_tot_p_nscf_rot(i,ik,3,1,:)),&	!fx,fy,fz summed over z
						(spin_kr_tot_p_nscf_rot(i,ik,1,0,0)), (spin_kr_tot_p_nscf_rot(i,ik,2,0,0)),(spin_kr_tot_p_nscf_rot(i,ik,3,0,0)),&	!cx,cy,cz z=0
						(spin_kr_tot_p_nscf_rot(i,ik,1,1,0)), (spin_kr_tot_p_nscf_rot(i,ik,2,1,0)),(spin_kr_tot_p_nscf_rot(i,ik,3,1,0))	!fx,fy,fz z=0
write(229,"(i5,3f10.5, f17.10,30e13.5)") ik, kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), kr_vecs_nscf(ik,3),e_kr_tot_nscf_pm(i,ik,2),&
						sum(spin_kr_tot_m_nscf_rot(i,ik,1,0,:)), sum(spin_kr_tot_m_nscf_rot(i,ik,2,0,:)),sum(spin_kr_tot_m_nscf_rot(i,ik,3,0,:)),&	!cx,cy,cz summed over z
						sum(spin_kr_tot_m_nscf_rot(i,ik,1,1,:)), sum(spin_kr_tot_m_nscf_rot(i,ik,2,1,:)),sum(spin_kr_tot_m_nscf_rot(i,ik,3,1,:)),&	!fx,fy,fz summed over z
						(spin_kr_tot_m_nscf_rot(i,ik,1,0,0)), (spin_kr_tot_m_nscf_rot(i,ik,2,0,0)),(spin_kr_tot_m_nscf_rot(i,ik,3,0,0)),&	!cx,cy,cz z=0
						(spin_kr_tot_m_nscf_rot(i,ik,1,1,0)), (spin_kr_tot_m_nscf_rot(i,ik,2,1,0)),(spin_kr_tot_m_nscf_rot(i,ik,3,1,0))	!fx,fy,fz z=0
						
     enddo
write(228,"(a10)") "   "
write(228,"(a10)") "   "
write(229,"(a10)") "   "
write(229,"(a10)") "   "
    enddo
 close(unit=227)






print *, "k dot p kr"
allocate(e_kr_cmplx(0:size_h_kr-1))
allocate(e_kr_tot_nscf_cmplx(0:size_h_kr-1,0:nkr_plot-1))
do ik=0, nkr_plot-1

!if (.not. l2dmodel .and. ik<=190) 			 call build_ham_kr_kdotp(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.)
!if (.not. l2dmodel .and. ik>190) 			 call build_ham_kr_kdotp_x(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.)

if (.not. l2dmodel) then
!theta,omega,phi,kx0,ky0
!if ((surf_type=="100" .or. surf_type=="001") .and. ik<=190)  call build_ham_kr_kdotp(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.)
if ((surf_type=="100" .or. surf_type=="001") .and. ik<=190) call build_ham_kr_kdotp_theta(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.,0.0d0,0.0d0,0.d0,0.d0,0.d0)	!G 001
if ((surf_type=="100" .or. surf_type=="001") .and. ik>190)  call build_ham_kr_kdotp_theta(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.,pi/2,0.d0,0.d0,0.5d0,0.d0)	!X 100
if ((surf_type=="111")) call build_ham_kr_kdotp_theta(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.,atan(sqrt(2.0d0)),0.0d0,pi/4,1./sqrt(6.d0),0.d0)			!M 
if ((surf_type=="110").and. ik<=170)call build_ham_kr_kdotp_theta(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.,pi/2,0.d0,0.0d0,0.5d0,0.0d0)				!X 110
if ((surf_type=="110").and. ik>170) call build_ham_kr_kdotp_theta(ik,Ham_kr,b_kr, lambda_kr,mu_kr,kr_vecs_nscf, nkr_plot, .true.,pi/4,pi/2,0.0d0,0.0d0,0.5d0/sqrt(2.))			!Y 101

endif


  call diagonalize_ham_kr(Ham_kr,e_kr,w_kr)
  e_kr_tot_nscf_cmplx(:,ik)=e_kr(:)
!  e_kr_tot_nscf(:,ik)=e_kr(:)
!  call diagonalize_ham_kr_kdotp(Ham_kr,e_kr_cmplx,w_kr)
!  e_kr_tot_nscf_cmplx(:,ik)=e_kr_cmplx(:)
enddo

deallocate(e_kr_cmplx)

OPEN(unit=227,file=trim(label)//'/bands_kr_kdotp',status='unknown')

    do i=0,size_h_kr-1
     do ik= 0, nkr_plot-1
write(227,"(i5,3f10.5, 2f20.10,15e12.5)") ik, kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), kr_vecs_nscf(ik,3),real(e_kr_tot_nscf_cmplx(i,ik)), aimag(e_kr_tot_nscf_cmplx(i,ik))!,&
!						(ew_kr_tot_nscf(i,ik,1,j),j=0,l_orb_tot-1),(ew_kr_tot_nscf(i,ik,2,j), j=0,l_orb_tot-1),(ew_kr_tot_nscf(i,ik,3,j), j=0,l_orb_tot-1)
     enddo
write(227,"(a10)") "   "
write(227,"(a10)") "   "
    enddo
 close (unit=227)


if (lpathnscf) then
call system('mkdir ' //trim(label)//"/wfc")
do ik=0, nkr_plot-1
 WRITE(UNIT=ak,FMT='(I4.4)') ik
OPEN(unit=228,file=trim(label)//'/wfc/wfc_kr_nscf_k'//ak,status='unknown')
do i=0,size_h_kr-1
do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
   ind1= index_xyzcso(ix,iy,iz,0,0,0,nx_kr,ny_kr,nz)   	!d1down
   ind2= index_xyzcso(ix,iy,iz,0,1,0,nx_kr,ny_kr,nz)   	!d1up
   ind3= index_xyzcso(ix,iy,iz,0,0,1,nx_kr,ny_kr,nz)   	!d2down
   ind4= index_xyzcso(ix,iy,iz,0,1,1,nx_kr,ny_kr,nz)   	!d2up
   ind5= index_xyzcso(ix,iy,iz,1,0,0,nx_kr,ny_kr,nz)   	!f1down
   ind6= index_xyzcso(ix,iy,iz,1,1,0,nx_kr,ny_kr,nz)   	!f1up
   ind7= index_xyzcso(ix,iy,iz,1,0,1,nx_kr,ny_kr,nz)   	!f2down
   ind8= index_xyzcso(ix,iy,iz,1,1,1,nx_kr,ny_kr,nz)   	!f2up
   ind9= index_xyzcso(ix,iy,iz,1,0,2,nx_kr,ny_kr,nz)   	!f3down
  ind10= index_xyzcso(ix,iy,iz,1,1,2,nx_kr,ny_kr,nz)   	!f3up

if (l_orbital==2 .and. l_orbital_f==2) then
write(228,"(i6,3i3,f12.5, 30e15.6)") , i,ix,iy,iz,e_kr_tot_nscf(i,ik),	&
			abs(w_kr_tot_nscf(ind1,i,ik))**2+abs(w_kr_tot_nscf(ind2,i,ik))**2,abs(w_kr_tot_nscf(ind3,i,ik))**2+abs(w_kr_tot_nscf(ind4,i,ik))**2,&
			abs(w_kr_tot_nscf(ind4,i,ik))**2+abs(w_kr_tot_nscf(ind5,i,ik))**2,abs(w_kr_tot_nscf(ind6,i,ik))**2+abs(w_kr_tot_nscf(ind7,i,ik))**2
endif


enddo
enddo
enddo
write(228,"(a1)"),""
enddo	!i eigenvalue
 close(unit=228)
 enddo


 endif

















deallocate(e_kr_tot_nscf_cmplx)
deallocate(ew_kr_tot_nscf)
deallocate(ew_kr_tot_nscf_rec)
deallocate(ew_kr_tot_p_nscf)
deallocate(ew_kr_tot_m_nscf)
deallocate(spin_kr_tot_nscf)
deallocate(spin_kr_tot_nscf_rot)
deallocate(spin_kr_tot_p_nscf)
deallocate(spin_kr_tot_m_nscf)
deallocate(spin_kr_tot_p_nscf_rot)
deallocate(spin_kr_tot_m_nscf_rot)

endif

!if (abs(T-T_min)<1e-6 .and. ldos) call compute_dos_kr_nscf



!deallocate(kr_vecs_nscf)
!deallocate(factor_kr_nscf)
!if  (.not. lreadham) deallocate(phi_kr_sigma_alpha_nscf)
!if  (.not. lreadham) deallocate(phi_kr_alpha_sigma_nscf)
!if  (.not. lreadham) deallocate(phi_kr_sigma_alpha2_nscf)
!if  (.not. lreadham) deallocate(phi_kr_alpha_sigma2_nscf)
!if  (.not. lreadham) deallocate(phi_kr_sigma_alpha3_nscf)
!if  (.not. lreadham) deallocate(phi_kr_alpha_sigma3_nscf)
!if  (.not. lreadham) deallocate(delta_kr_nscf)
!if  (.not. lreadham) deallocate(kin_energy_kr_nscf)
!if  (.not. lreadham) deallocate(kin_energy2_kr_nscf)
!if  (.not. lreadham) deallocate(kin_energy3_kr_nscf)
!deallocate(corr_fct_kr_nscf)
!deallocate(v_kr_nscf)
!deallocate(fermien_kr)
!deallocate(occup_kr)
if (.not. ldiag_g0k) deallocate(e_kr_tot_nscf)
if (.not. ldiag_g0k) deallocate(w_kr_tot_nscf)
!if (ldiag_g0k) deallocate(H_kr_tot_nscf)

if (allocated(m_xy)) deallocate(m_xy)







end subroutine do_kr_nscf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kvecs (k_vecs)
!real(kind=idp),intent(out)	::k_vecs(0:n_k-1,1:4)
real(kind=idp),intent(out),allocatable	::k_vecs(:,:)

integer	::ikx,iky,ikz

allocate(k_vecs(0:n_k-1,1:4))

k_vecs=0

!simple cubic

do ik=0,n_k-1
   
   call convert_index_k(ik,ikx,iky,ikz)
   
if ((mod(nkx_k,2)==0 .and. .not. lkshifted) .or. (mod(nkx_k,2)==1 .and.  lkshifted))   then
k_vecs (ik,1)=-0.5+real(ikx)/real(nkx_k)
else
k_vecs (ik,1)=-0.5+real(ikx+0.5)/real(nkx_k)	!k_x 	odd
endif


if ((mod(nky_k,2)==0 .and. .not. lkshifted) .or. (mod(nky_k,2)==1 .and.  lkshifted))   then
k_vecs (ik,2)=-0.5+real(iky)/real(nky_k)
else
k_vecs (ik,2)=-0.5+real(iky+0.5)/real(nky_k)	!k_x 	odd
endif


!z: sempre non shifted
if (mod(nkz,2)==1)   k_vecs (ik,3)=-0.5+real(ikz+0.5)/real(nkz)	!k_z
if (mod(nkz,2)==0)   k_vecs (ik,3)=-0.5+real(ikz)/real(nkz)	!k_z

  
enddo

k_vecs(:,4)=1.0d0



end subroutine set_kvecs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_kvecs2d (k_vecs,nkvecs, nkx, nky, igamma)
!real(kind=idp),intent(out)	::k_vecs(0:n_k-1,1:4)
real(kind=idp),intent(out),allocatable	::k_vecs(:,:)
integer, intent(in)	::nkx, nky
integer	::ikx,iky,ikz
integer,intent(out)	::nkvecs
logical	::igamma


nkvecs=nkx*nky

allocate(k_vecs(0:nkvecs-1,1:4))

k_vecs=0

!simple cubic: gamma must be in

do ikx=0,nkx-1
 do iky=0,nky-1
   ik=ikx+nkx*iky
   
if (igamma .and. mod(nkx,2)==1)   k_vecs (ik,1)=-0.5d0+real(ikx+0.5d0)/real(nkx)	!k_x 	odd
if (igamma .and. mod(nky,2)==1)   k_vecs (ik,2)=-0.5d0+real(iky+0.5d0)/real(nky)	!k_y

if (igamma .and. mod(nkx,2)==0)   k_vecs (ik,1)=-0.5d0+real(ikx)/real(nkx)	!k_x	even
if (igamma .and. mod(nky,2)==0)   k_vecs (ik,2)=-0.5d0+real(iky)/real(nky)	!k_y


if (.not. igamma .and. mod(nkx,2)==0)   k_vecs (ik,1)=-0.5d0+real(ikx+0.5d0)/real(nkx)	!k_x 	odd
if (.not. igamma .and. mod(nky,2)==0)   k_vecs (ik,2)=-0.5d0+real(iky+0.5d0)/real(nky)	!k_y

if (.not. igamma .and. mod(nkx,2)==1)   k_vecs (ik,1)=-0.5d0+real(ikx)/real(nkx)	!k_x	even
if (.not. igamma .and. mod(nky,2)==1)   k_vecs (ik,2)=-0.5d0+real(iky)/real(nky)	!k_y

  enddo
enddo

k_vecs(:,3)=0.0d0	!z component
k_vecs(:,4)=1.0d0	!weight

end subroutine set_kvecs2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kvecs2d_tr (k_vecs,nkvecs, nkx, nky, igamma)
!real(kind=idp),intent(out)	::k_vecs(0:n_k-1,1:4)
real(kind=idp),intent(out),allocatable	::k_vecs(:,:)
integer, intent(in)	::nkx, nky
integer	::ikx,iky,ikz
integer,intent(out)	::nkvecs
logical	::igamma
integer	::nky_new

!takes into account time reversal
if (igamma) 		nky_new=int((nky+2)/2)
if (.not. igamma) 	nky_new=int((nky+1)/2)

nkvecs=nkx*nky_new

allocate(k_vecs(0:nkvecs-1,1:4))

k_vecs=0

!simple cubic: gamma must be in

do ikx=0,nkx-1
 do iky=0,nky_new-1
   ik=ikx+nkx*iky
   
if (igamma .and. mod(nkx,2)==1)   k_vecs (ik,1)=-0.5+real(ikx+0.5)/real(nkx)	!k_x 	odd
if (igamma .and. mod(nky,2)==1)   k_vecs (ik,2)=real(iky)/real(nky)	!k_y
if (igamma .and. mod(nky,2)==1 .and. iky==0) k_vecs(ik,4)=1.0d0
if (igamma .and. mod(nky,2)==1 .and. iky>0) k_vecs(ik,4)=2.0d0


if (igamma .and. mod(nkx,2)==0)   k_vecs (ik,1)=-0.5+real(ikx)/real(nkx)	!k_x	even
if (igamma .and. mod(nky,2)==0)   k_vecs (ik,2)=real(iky)/real(nky)	!k_y
if (igamma .and. mod(nky,2)==0 .and. iky==0) k_vecs(ik,4)=1.0d0
if (igamma .and. mod(nky,2)==0 .and. iky>0) k_vecs(ik,4)=2.0d0
if (igamma .and. mod(nky,2)==0 .and. iky==nky_new-1) k_vecs(ik,4)=1.0d0


if (.not. igamma .and. mod(nkx,2)==0)   k_vecs (ik,1)=-0.5+real(ikx+0.5)/real(nkx)	!k_x 	odd
if (.not. igamma .and. mod(nky,2)==0)   k_vecs (ik,2)=real(iky+0.5)/real(nky)	!k_y
if (.not. igamma .and. mod(nky,2)==0) k_vecs(ik,4)=2.0d0



if (.not. igamma .and. mod(nkx,2)==1)   k_vecs (ik,1)=-0.5+real(ikx)/real(nkx)	!k_x	even
if (.not. igamma .and. mod(nky,2)==1)   k_vecs (ik,2)=real(iky+0.5)/real(nky)	!k_y
if (.not. igamma .and. mod(nky,2)==1 .and. iky<nky_new-1) k_vecs(ik,4)=2.0d0
if (.not. igamma .and. mod(nky,2)==1 .and. iky==nky_new-1) k_vecs(ik,4)=1.0d0



  enddo
enddo

k_vecs(:,3)=0.0d0	!z component
!k_vecs(:,4)=1.0d0	!weight

end subroutine set_kvecs2d_tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_v_kr(v_kr_tot)
complex(kind=idp),intent(out):: v_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1)
integer	::ik

print *, "build_v_kr"

v_kr_tot=0

do ik=0, nk2d-1
 
!if (.not. lreadham)  call build_ham_kr(ik,V_kr,one_kr,zero_kr, 0.0d0,.false.,nk2d,phi_kr_sigma_alpha_tot,phi_kr_alpha_sigma_tot, kin_energy_kr_tot, kin_energy2_kr_tot, kin_energy3_kr_tot,&	
!  			phi_kr_sigma_alpha2_tot,phi_kr_alpha_sigma2_tot,phi_kr_sigma_alpha3_tot,phi_kr_alpha_sigma3_tot,k_vecs2d)
!if (lreadham)        call build_ham_kr_read(ik,V_kr,one_kr,zero_kr, 0.0d0,.false.,nk2d,k_vecs2d)

call build_ham_kr2(ik,V_kr,one_kr, zero_kr,0.0d0,k_vecs2d,nk2d, .false.)

v_kr_tot(:,:,ik)=v_kr(:,:)
  
enddo

end subroutine build_v_kr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_V_kr2(h_kr_psi_tot, index_kr_psi_tot)
integer			::index_kr_psi_tot(0:size_h_kr-1,n_max_link,3)
complex(kind=idpc),allocatable	::h_kr_psi_tot(:,:,:)

print *, "build_v_kr2"

allocate(V(0:size_h_kr-1,0:size_h_kr-1))
allocate(h_kr_psi_tot(0:size_h_kr-1,n_max_link,0:nk2d-1))
!allocate(w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1))
!allocate(e_kr_tot(0:size_h_kr-1,0:nk2d-1))


do ik=0, nk2d-1

 V=0
 
 do ix=0,nx_kr-1	!x
  do iy=0,ny_kr-1	!y
   do iz=0,nz-1	!z
    do icf=0,1	!c/f
     do is=0,1	!up/down 
      do o1=0,l_orb(icf)-1
      
      indp=index_xyzcso (ix , iy , iz, icf, is, o1, nx_kr, ny_kr, nz)
      
      call ham_psi(ik,ix,iy,iz,icf,is,o1,index_psi,h_psi,one_kr,zero_kr,0.0d0,.false.,nx_kr,ny_kr,nz,ec_site_kr,ef_site_kr, k_vecs2d,nk2d)

      do i=1,n_max_link
       ind2p=index_psi(i,7)
       V(indp,ind2p)=V(indp,ind2p)+h_psi(i)
       h_kr_psi_tot(indp,i,ik)=h_psi(i)
       index_kr_psi_tot(indp,i,1)=index_psi(i,6)	!ind2	xyz
       index_kr_psi_tot(indp,i,2)=index_psi(i,7)	!ind2p	xyzcso
       index_kr_psi_tot(indp,i,3)=index_psi(i,4)	!icf	a
      enddo
      
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

 do i=0, size_h_kr-1
 do j=i, size_h_kr-1
 if (abs(V(j,i)-conjg(V(i,j)))>1e-9) then		! v is only for check reasons
 print *, "error in V", ik,i, j, V(i,j), V(j,i)
 !stop
 endif
 enddo
 enddo

enddo	!k


deallocate(v)

end subroutine build_V_kr2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_initialvalues_k

INQUIRE(file=trim(label)//'/b_k2', EXIST=lfileexists)
If(lfileexists) print *, "File b_k2 found"
If(.not. lfileexists) print *, "File b_k2 not found"


if (lreadbk .and. lfileexists) then
 write (*, "(a40)") "Reading b, lambda, mu from file"
 OPEN(unit=30,file=trim(label)//'/b_k2',status='unknown')
 read(30,"(3f20.10)"), b_k, lambda_k, mu_k		!file b_k2
 close (unit=30)
else
 b_k=b_init
 lambda_k=lambda_init
 mu_k=mu_init
endif

write (*, "(a30)") "Initial b, lambda, mu"
write(*,"(3f20.10)"), b_k, lambda_k, mu_k

!on site energies
if (l_orbital>0) ec_site_k(0)=ec_
if (l_orbital>1) ec_site_k(1)=ec_2
if (l_orbital_f>0) ef_site_k(0)=ef_
if (l_orbital_f>1) ef_site_k(1)=ef_2
if (l_orbital_f>2) ef_site_k(2)=ef_7



end subroutine set_initialvalues_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_initialvalues_kr
integer		::nz2, nx2, ny2
real(kind=idp)	::mu_temp

print *, "set_initialvalues_kr"

lvac_kr=0

INQUIRE(file=trim(label)//'/b_kr2', EXIST=lfileexists)

if (loopk) then				!if not, uses values from k loop
 b_kr=b_k
 lambda_kr=lambda_k
 mu_kr=mu_k
endif

if (lreadbkr .and. lfileexists) then				!reads from file
 write (*, "(a40)") "Reading b_kr, lambda_kr, mu_kr from file"

 OPEN(unit=31,file=trim(label)//'/b_kr2',status='unknown')

     read(31,"(3i4,100f15.10)") nx2, ny2, nz2
 do ix=0,nx2-1	
   do iy=0,ny2-1
    do iz=0,nz2-1

!     ind=ix + iy*nx + iz*nx*ny		
     read(31,"(3i4,100f15.10)") ix2,iy2,iz2, b_kr(indxx2(ix2,nx2,iy2,ny2,iz2,nz2, nx_kr, ny_kr) ), lambda_kr(indxx2(ix2,nx2,iy2,ny2,iz2,nz2,nx_kr,ny_kr))
 

   enddo
  enddo
 enddo
 
 read(31,"(f20.10)"), mu_temp
 
 write (*, "(a40)") "OK Reading b_kr, lambda_kr, mu_kr"


 if (.not. lfixmu_kr) mu_kr=mu_temp
 close (unit=31)

 do ix=0,nx_kr-1	
   do iy=0,ny_kr-1
    do iz=0,nz-1

!     ind=ix + iy*nx + iz*nx*ny		
ind=index_xyz(ix,iy,iz,nx_kr,ny_kr) 
ind2=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr) 
!print *, ix,iy,iz,ind,ind2
b_kr(ind)=b_kr(ind2)
lambda_kr(ind)=lambda_kr(ind2)

   enddo
  enddo
 enddo


! if(nz==nz2) mu_kr=mu_temp			!reads mu only if I have the same system size, otherwise use value from k loop

else if (loopk) then				!if not, uses values from k loop
 b_kr=b_k
 lambda_kr=lambda_k
 mu_kr=mu_k
else						!if not, using starting values
 b_kr=b_init
 lambda_kr=lambda_init
 mu_kr=mu_init
endif



if (l_orbital_f>0) ef_site_kr(:,0)=ef_
if (l_orbital_f>1) ef_site_kr(:,1)=ef_2
if (l_orbital_f>2) ef_site_kr(:,2)=ef_7
if (l_orbital>0) ec_site_kr(:,0)=ec_
if (l_orbital>1) ec_site_kr(:,1)=ec_2

!energies
do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1
 indp=index_xyz (ix , iy , iz,nx_kr, ny_kr)


!if (l_orbital_f>0) ef_site_kr(indp,0)=ef_site_kr(indp,0)+d_ef_kr*(exp(-iz/lambda_dec_kr)+exp((iz-nz+1)/lambda_dec_kr))
!if (l_orbital_f>1) ef_site_kr(indp,1)=ef_site_kr(indp,1)+d_ef_kr*(exp(-iz/lambda_dec_kr)+exp((iz-nz+1)/lambda_dec_kr))
!if (l_orbital_f>2) ef_site_kr(indp,2)=ef_site_kr(indp,2)+d_ef7_kr*(exp(-iz/lambda_dec_kr)+exp((iz-nz+1)/lambda_dec_kr))
!if (l_orbital>0)   ec_site_kr(indp,:)=ec_site_kr(indp,:)+d_ec_kr*(exp(-iz/lambda_dec_kr)+exp((iz-nz+1)/lambda_dec_kr))


!scattering pot on first and last layer
if(iz==0 .or. iz==nz-1) then
if (l_orbital_f>0) ef_site_kr(indp,0)=ef_site_kr(indp,0)+d_ef_kr
if (l_orbital_f>1) ef_site_kr(indp,1)=ef_site_kr(indp,1)+d_ef_kr_82
if (l_orbital_f>2) ef_site_kr(indp,2)=ef_site_kr(indp,2)+d_ef7_kr
if (l_orbital>0)   ec_site_kr(indp,:)=ec_site_kr(indp,:)+d_ec_kr
endif

if(nz>2 .and. (iz==1 .or. iz==nz-2)) then
if (l_orbital_f>0) ef_site_kr(indp,0)=ef_site_kr(indp,0)+d_ef_2_kr
if (l_orbital_f>1) ef_site_kr(indp,1)=ef_site_kr(indp,1)+d_ef_2_kr
if (l_orbital_f>2) ef_site_kr(indp,2)=ef_site_kr(indp,2)+d_ef7_2_kr
if (l_orbital>0)   ec_site_kr(indp,:)=ec_site_kr(indp,:)+d_ec_2_kr
endif


!remove an atom for reconstruction
if((iz==0 .or. iz==nz-1) .and. nx_kr==2 .and. ix==1 .and. lrec) then
if (l_orbital_f>0) ef_site_kr(indp,0)=ef_site_kr(indp,0)+d_e_rec			!+1e2
if (l_orbital_f>1) ef_site_kr(indp,1)=ef_site_kr(indp,1)+d_e_rec	!+1e2
if (l_orbital_f>2) ef_site_kr(indp,2)=ef_site_kr(indp,2)+d_e_rec	!+1e2
if (l_orbital>0)   ec_site_kr(indp,:)=ec_site_kr(indp,:)+d_e_rec	!+1e2
lvac_kr(indp)=1	!ok for scf loop
print *, "rec",indp!,lvac_kr(indp)
endif

if((iz==1 .or. iz==nz-2) .and. nx_kr==2 .and. ix==1 .and. lrec) then
if (l_orbital_f>0) ef_site_kr(indp,0)=ef_site_kr(indp,0)+d_e_rec2			!+1e2
if (l_orbital_f>1) ef_site_kr(indp,1)=ef_site_kr(indp,1)+d_e_rec2	!+1e2
if (l_orbital_f>2) ef_site_kr(indp,2)=ef_site_kr(indp,2)+d_e_rec2	!+1e2
if (l_orbital>0)   ec_site_kr(indp,:)=ec_site_kr(indp,:)+d_e_rec2	!+1e2
!lvac_kr(indp)=1	!ok for scf loop
print *, "rec2",indp!,lvac_kr(indp)
endif



!ef_site_kr(nz-1,0)=ef_site_kr(nz-1,0)+d_ef_kr
!if (l_orbital_f>1) ef_site_kr(nz-1,1)=ef_site_kr(nz-1,1)+d_ef_kr
!if (l_orbital_f>2) ef_site_kr(nz-1,2)=ef_site_kr(nz-1,2)+d_ef7_kr
!ec_site_kr(nz-1,:)=ec_site_kr(nz-1,:)+d_ec_kr

!if (limpurity_kr) ef_site_kr(:,0)=ef_site_kr(:,0)+ef0
!if (limpurity_kr) ec_site_kr(0)=ec_site_kr(0)+d_ec_kr
!if (limpurity_kr) ef_site_kr(:,nz-1)=ef_site_kr(:,nz-1)+ef0
!if (limpurity_kr) ec_site_kr(nz-1)=ec_site_kr(nz-1)+ec0
!if (limpurity_kr2) ef_site_kr(nz-1)=ef_site_kr(:,nz-1)
!if (limpurity_kr2) ec_site_kr(nz-1)=ec0

enddo
enddo
enddo


!vacancies

do ix=0, nx_kr-1
do iy=0, ny_kr-1
do iz=0, nz-1

 indp=index_xyz (ix , iy , iz,nx, ny)

if (lvacancy_kr .and. ix==0 .and. iy==0 .and. (iz==0 .or. iz==nz-1)) then
 ef_site_kr(indp,:)=1e2
 ec_site_kr(indp,:)=1e2
 lvac_kr(indp)=1
endif

enddo
enddo
enddo


!if (lvacancy_kr2) ef_site_kr(nz-1,:)=1e2
!if (lvacancy_kr3) ef_site_kr(nz/2,:)=1e2
!if (lvacancy_kr4) ef_site_kr(1,:)=1e2
!if (lvacancy_kr5) ef_site_kr(nz-2,:)=1e2
!if (lvacancy_kr6) ef_site_kr(nz/2-1,:)=1e2

!if (lvacancy_kr)  lvac_kr(0)=1
!if (lvacancy_kr2) lvac_kr(nz-1)=1
!if (lvacancy_kr3) lvac_kr(nz/2)=1
!if (lvacancy_kr4) lvac_kr(1)=1
!if (lvacancy_kr5) lvac_kr(nz-2)=1
!if (lvacancy_kr6) lvac_kr(nz/2-1)=1

do i=0,n_sites_kr-1
!if (lvac_kr(i)) b_kr(i)=0
!if (lvac_kr(i)) lambda_kr(i)=100
enddo
!b_kr(0)=0.02
!b_kr(nz-1)=0.02

!n el
print *, "n_el_kr without holes", n_el_kr 
do i=0, nz-1
if(lvac_kr(i)) N_el_kr=N_el_kr-1
enddo
print *, "n_el_kr with holes", n_el_kr 




!OPEN(unit=100,file=trim(label)//'/eb_kr_in',status='unknown')

write (*, "(a40)") "Initial b_kr, lambda_kr, ef, ec, vac"!, vac, mu_kr:"
do ix=0,nx_kr-1
 do iy=0,ny_kr-1
  do iz=0,nz-1
   i=index_xyz(ix,iy,iz,nx_kr,ny_kr)
 write(*,"(3i5,20f17.10)"), ix,iy,iz,b_kr(i), lambda_kr(i), (ef_site_kr(i,j), j=0,l_orbital_f-1),  (ec_site_kr(i,j), j=0,l_orbital-1), dble(lvac_kr(i))
! write(100,"(4f20.10)"), b_kr(i), lambda_kr(i), ef_site_kr(i), ec_site_kr(i)
enddo
enddo
enddo
write(*,"(f20.10)"), mu_kr

! close(unit=100) 


end subroutine set_initialvalues_kr

!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_initialvalues_r
integer	:: nx2, ny2, nz2, index_kr
real(kind=idp)	:: mu_temp

if(l2dmodel) then
b=1
lambda=0
mu=0
ef_site=0
ec_site=0
lvac_r=0
endif

if (.not. l2dmodel) then

INQUIRE(file=trim(label)//'/b_r2', EXIST=lfileexists)

if (loopkr) then
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
    ind=ix + iy*nx + iz*nx*ny		
    b(ind)=b_new_kr(iz)
    lambda(ind)=lambda_new_kr(iz)
   enddo
  enddo
 enddo
 mu=mu_new_kr
endif


if (lreadbr .and. lfileexists) then				!reads from file
 write (*, "(a40)") "Reading b_r, lambda_r, mu_r from file"

 OPEN(unit=37,file=trim(label)//'/b_r2',status='unknown')

     read(37,"(3i4,100f15.10)") nx2, ny2, nz2
 do ix=0,nx2-1	
   do iy=0,ny2-1
    do iz=0,nz2-1

!     ind=ix + iy*nx + iz*nx*ny		
     read(37,"(3i4,100f15.10)") ix2,iy2,iz2, b(indxx2(ix2,nx2,iy2,ny2,iz2,nz2,nx,ny) ), lambda(indxx2(ix2,nx2,iy2,ny2,iz2,nz2,nx,ny))
 

   enddo
  enddo
 enddo

 read(37,"(f20.10)")  mu_temp
 
 if (nx2==nx .and. ny2==ny .and. nz2==nz) mu=mu_temp	!reads mu only if same system size
 if (loopkr) mu=mu_new_kr

 close (unit=37)

else if (loopkr) then
 write (*, "(a40)") "Taking b_r, lambda_r, mu_r from kr loop"
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
    ind=index_xyz(ix,iy,iz,nx,ny)
    index_kr=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr)
    b(ind)=b_kr(index_kr)
    lambda(ind)=lambda_kr(index_kr)
   enddo
  enddo
 enddo
 mu=mu_new_kr
else if (loopk) then
 write (*, "(a40)") "Taking b_r, lambda_r, mu_r from k loop"
 b=b_new_k
 lambda=lambda_new_k
 mu=mu_new_k
else
 write (*, "(a40)") "Taking b_r, lambda_r, mu_r from input"
 b=b_init
 lambda=lambda_init
 mu=mu_init
endif

!vacancies
do i=0, n_sites-1
if (lvac_r(i)) b(i)=0
if (lvac_r(i)) lambda(i)=-10
enddo


print *, nx,ny,nz


     write(*,"(a42)") "Initial b,lambda, ef,ec,vac"
 do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z

     ind=ix + iy*nx + iz*nx*ny		
     write(*,"(3i4,10f15.10)") ix,iy,iz, b(ind) , lambda(ind),  (ef_site(ind,j), j=0,l_orbital_f-1),(ec_site(ind,j), j=0,l_orbital-1),dble(lvac_r(ind))
 

   enddo
  enddo
 enddo
 write(*,"(a10,f20.10)") "mu: ", mu

endif !l2dmodel


!stop

end subroutine set_initialvalues_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_siteenergies_r
real(kind=idp)	::delta
integer	::nx2,ny2,nz2

if(loopkr) then
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1


ind=index_xyz(ix,iy,iz,nx,ny)
ind2=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr)
ef_site(ind,:)=ef_site_kr(ind2,:)
ec_site(ind,:)=ec_site_kr(ind2,:)

enddo
enddo
enddo
endif

INQUIRE(file=trim(label)//'/ef_r2', EXIST=lfileexists)

!if (lreadbr .and. lfileexists) then				!reads from file

!lvac_r=0

! write (*, "(a40)") "Reading ef_r, ec from file"

! OPEN(unit=37,file=trim(label)//'/ef_r2',status='unknown')

!     read(37,"(3i4,100f15.10)") nx2, ny2, nz2
     
! do ix=0,nx2-1	
!   do iy=0,ny2-1
!    do iz=0,nz2-1

!     ind=ix + iy*nx + iz*nx*ny		
!if (.not. lmagimp) then
! read(37,"(3i4,100e15.6)") ix2,iy2,iz2, ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0),ec_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0)
 !ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),1)=ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0)
 !ec_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),1)=ec_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0)
!else
! read(37,"(3i4,100e15.6)") ix2,iy2,iz2, (ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),j),j=0,l_orbital_f-1),&
! ec_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0)!,ec_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),1)
!endif
!     if (ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2),0)>0.99*1e2) lvac_r(indxx2(ix2,nx2,iy2,ny2,iz2,nz2))=1	!finds vacancies
 
!      write(*,"(3i4,100e15.6)") ix2,iy2,iz2, ef_site(indxx2(ix2,nx2,iy2,ny2,iz2,nz2) )

!   enddo
!  enddo
! enddo


! close (unit=37)



!else

lvac_r=0

!if (mod(nx,2)==1) ind_vac=(nx/2) + (ny/2)*nx
!if (mod(nx,2)==1) ind_vac2=index_xyz(nx/2,ny/2,0,nx,ny)

!if (mod(nx,2)==0) ind_vac=(nx/2)-1 + (ny/2)*nx
!if (mod(nx,2)==0) ind_vac2=index_xyz(nx/2-1,ny/2,0,nx,ny)
!ind_vac=0
ind_vac =index_xyz((nx-1)/2,(ny-1)/2,0,nx,ny)
!ind_vac2=(nx/2) + (ny/2)*nx + nx*ny *(0)
ind_vac2=index_xyz((nx-1)/2+1,(ny-1)/2,0,nx,ny)
!ind_vac2=(nx/2) + (ny/2)*nx + nx*ny *(nz-1)
!ind_vac3=(nx/2) + (ny/2)*nx + nx*ny *(nz/2)
ind_vac3=index_xyz((nx-1)/2,(ny-1)/2,nz/2,nx,ny)
!ind_vac4=(nx/2) + (ny/2)*nx + nx*ny *(1)
ind_vac4=index_xyz((nx-1)/2,(ny-1)/2,1,nx,ny)
!ind_vac5=(nx/2) + (ny/2)*nx + nx*ny *(nz-2)
ind_vac5=index_xyz((nx-1)/2,(ny-1)/2,nz-2,nx,ny)

!ind=ix + iy*nx + iz*nx*ny		


 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1

ind=index_xyz(ix,iy,iz,nx,ny)		
			        
!if (.not. lmagimp) then
!1st layer
!stripes of impurities
!if (limpurity .and. mod(ix,nx_kr)==0 .and. iz==0) ef_site(ind,0)=ef_site(ind,0)+d_ef
!if (limpurity .and. mod(ix,nx_kr)==0 .and. iz==0 .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
!if (limpurity .and. mod(ix,nx_kr)==0 .and. iz==0 .and. l_orbital_f>2 ) ef_site(ind,2)=ef_site(ind,2)+d_ef7
!if (limpurity .and. mod(ix,nx_kr)==0 .and. iz==0 ) ec_site(ind,:)=ec_site(ind,:)+d_ec

!if (limpurity2 .and. mod(ix,nx_kr)==0 .and. iz==1) ef_site(ind,0)=ef_site(ind,0)+d_ef
!if (limpurity2 .and. mod(ix,nx_kr)==0 .and. iz==1 .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
!if (limpurity2 .and. mod(ix,nx_kr)==0 .and. iz==1 .and. l_orbital_f>2 ) ef_site(ind,2)=ef_site(ind,2)+d_ef7
!if (limpurity2 .and. mod(ix,nx_kr)==0 .and. iz==1 ) ec_site(ind,:)=ec_site(ind,:)+d_ec

if (limpurity .and. ind==ind_vac ) ef_site(ind,0)=ef_site(ind,0)+d_ef
if (limpurity .and. ind==ind_vac .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
if (limpurity .and. ind==ind_vac .and. l_orbital_f>2 ) ef_site(ind,2)=ef_site(ind,2)+d_ef7
if (limpurity .and. ind==ind_vac ) ec_site(ind,:)=ec_site(ind,:)+d_ec
!2nd layer
if (limpurity2 .and. ind==ind_vac4 ) ef_site(ind,0)=ef_site(ind,0)+d_ef
if (limpurity2 .and. ind==ind_vac4  .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
if (limpurity2 .and. ind==ind_vac4  .and. l_orbital_f>2) ef_site(ind,2)=ef_site(ind,2)+d_ef7
if (limpurity2 .and. ind==ind_vac4 ) ec_site(ind,:)=ec_site(ind,:)+d_ec
!last layer
!first layer x==1
!stripes of impurities
!if (limpurity2 .and. mod(ix,nx_kr)==1 .and. iz==0 ) ef_site(ind,0)=ef_site(ind,0)+d_ef
!if (limpurity2 .and. mod(ix,nx_kr)==1 .and. iz==0 .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
!if (limpurity2 .and. mod(ix,nx_kr)==1 .and. iz==0 .and. l_orbital_f>2 ) ef_site(ind,2)=ef_site(ind,2)+d_ef7
!if (limpurity2 .and. mod(ix,nx_kr)==1 .and. iz==0 ) ec_site(ind,:)=ec_site(ind,:)+d_ec


!if (limpurity2 .and. ind==ind_vac2 ) ef_site(ind,0)=ef_site(ind,0)+d_ef
!if (limpurity2 .and. ind==ind_vac2  .and. l_orbital_f>1) ef_site(ind,1)=ef_site(ind,1)+d_ef
!if (limpurity2 .and. ind==ind_vac2  .and. l_orbital_f>2) ef_site(ind,2)=ef_site(ind,2)+d_ef7
!if (limpurity2 .and. ind==ind_vac2 ) ec_site(ind,:)=ec_site(ind,:)+d_ec





if (lvacancy .and. ind==ind_vac )    ef_site(ind,:)=ef_site(ind,:)+1e2	!both Gamma7 and Gamma8
if (lvacancy2 .and. ind==ind_vac2 )  ef_site(ind,:)=ef_site(ind,:)+1e2
if (lvacancy3 .and. ind==ind_vac3 )  ef_site(ind,:)=ef_site(ind,:)+1e2
if (lvacancy4 .and. ind==ind_vac4 )  ef_site(ind,:)=ef_site(ind,:)+1e2
if (lvacancy5 .and. ind==ind_vac5 )  ef_site(ind,:)=ef_site(ind,:)+1e2

!if (lvacancy .and. ind==ind_vac )    ec_site(ind,:)=ec0
!if (lvacancy2 .and. ind==ind_vac2 )  ec_site(ind,:)=ec0
!if (lvacancy3 .and. ind==ind_vac3 )  ec_site(ind,:)=ec0
!if (lvacancy4 .and. ind==ind_vac4 )  ec_site(ind,:)=ec0
!if (lvacancy5 .and. ind==ind_vac5 )  ec_site(ind,:)=ec0


if (lvacancy .and. ind==ind_vac )   lvac_r(ind)=1
if (lvacancy2 .and. ind==ind_vac2 ) lvac_r(ind)=1
if (lvacancy3 .and. ind==ind_vac3 ) lvac_r(ind)=1
if (lvacancy4 .and. ind==ind_vac4 ) lvac_r(ind)=1
if (lvacancy5 .and. ind==ind_vac5 ) lvac_r(ind)=1
if (limpurity .and. ind==ind_vac .and. abs(d_ef+ef_)>50  .and.  abs(d_ef7+ef_7)>50)  lvac_r(ind)=1
if (limpurity2 .and. ind==ind_vac2 .and. abs(d_ef+ef_)>50  .and.  abs(d_ef7+ef_7)>50)  lvac_r(ind)=1


!else

!if (limpurity .and. ind==ind_vac ) ef_site(ind,0)=ef0_up
!if (limpurity .and. ind==ind_vac ) ec_site(ind,0)=ec0_up
!if (limpurity2 .and. ind==ind_vac4 ) ef_site(ind,0)=ef0_up
!if (limpurity2 .and. ind==ind_vac4 ) ec_site(ind,0)=ec0_up

!if (lvacancy .and. ind==ind_vac )    ef_site(ind,0)=1e2
!if (lvacancy2 .and. ind==ind_vac2 )  ef_site(ind,0)=1e2
!if (lvacancy3 .and. ind==ind_vac3 )  ef_site(ind,0)=1e2
!if (lvacancy4 .and. ind==ind_vac4 )  ef_site(ind,0)=1e2
!if (lvacancy5 .and. ind==ind_vac5 )  ef_site(ind,0)=1e2


!if (limpurity .and. ind==ind_vac ) ef_site(ind,1)=ef0_down
!if (limpurity .and. ind==ind_vac ) ec_site(ind,1)=ec0_down
!if (limpurity2 .and. ind==ind_vac4 ) ef_site(ind,1)=ef0_down
!if (limpurity2 .and. ind==ind_vac4 ) ec_site(ind,1)=ec0_down



!if (lvacancy .and. ind==ind_vac )    ef_site(ind,1)=-1e2
!if (lvacancy2 .and. ind==ind_vac2 )  ef_site(ind,1)=-1e2
!if (lvacancy3 .and. ind==ind_vac3 )  ef_site(ind,1)=-1e2
!if (lvacancy4 .and. ind==ind_vac4 )  ef_site(ind,1)=-1e2
!if (lvacancy5 .and. ind==ind_vac5 )  ef_site(ind,1)=-1e2


!if (lvacancy .and. ind==ind_vac )   lvac_r(ind)=1
!if (lvacancy2 .and. ind==ind_vac2 ) lvac_r(ind)=1
!if (lvacancy3 .and. ind==ind_vac3 ) lvac_r(ind)=1
!if (lvacancy4 .and. ind==ind_vac4 ) lvac_r(ind)=1
!if (lvacancy5 .and. ind==ind_vac5 ) lvac_r(ind)=1

!endif


   enddo
  enddo
 enddo


print *, lvac_r
    call random_seed

if(ldisorder_r) then
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1, max(nz-1,1)		!first and last layer

    call random_number(delta)


 ! print *, ix,iy,iz,delta
 
    ind=ix + iy*nx + iz*nx*ny		

   if(delta<conc_vac) then
    ef_site(ind,:)=ef_site(ind,:)+1e2
    lvac_r(ind)=1
   endif

   enddo
  enddo
 enddo
endif


if(ldisorder_r2) then
 do ix=0,nx-1
  do iy=0,ny-1

    call random_number(delta)

    iz=nz/2		!middle layer

 ! print *, ix,iy,iz,delta
 
    ind=ix + iy*nx + iz*nx*ny		

   if(delta<conc_vac) then
    ef_site(ind,:)=1e2
    lvac_r(ind)=1
   endif

  enddo
 enddo
endif


if(ldisorder_r_bulk) then
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=1,nz-2		!all layers except surfaces
    call random_number(delta)

 
    ind=ix + iy*nx + iz*nx*ny		

   if(delta<conc_vac) then
    ef_site(ind,:)=1e2
    lvac_r(ind)=1
   endif

  enddo
 enddo
enddo
endif






!writes to file

 OPEN(unit=37,file=trim(label)//'/ef_r2',status='unknown')

 write(37,"(3i4)") nx,ny,nz

 do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z

 ind=ix + iy*nx + iz*nx*ny		
!if (.not. lmagimp) write(37,"(3i4,100e15.6)") ix,iy,iz, ef_site(ind,0), ec_site(ind,0)
!if (lmagimp) 
write(37,"(3i4,100e15.6)") ix,iy,iz, (ef_site(ind,j), j=0,l_orbital_f-1) ,(ec_site(ind,j), j=0,l_orbital-1)!,ec_site(ind,1)

   enddo
  enddo
 enddo

 close (unit=37)

	
!endif


!stop


! sets number of electrons

print *, "N_el witohout holes=", N_el


 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1


ind=ix + iy*nx + iz*nx*ny		
if (lvac_r(ind)) N_el=N_el-1

enddo
enddo
enddo

print *, "N_el with holes=", N_el



end subroutine set_siteenergies_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function indxx2(ix2,nx2,iy2,ny2,iz2,nz2,nx,ny)
integer, intent(in)	::ix2,nx2,iy2,ny2,iz2,nz2,nx,ny
integer	::indx, indy, indz


!if (ix2< (nx2-1)/2) then
!indx=ix2
!else
!indx=ix2+nx-nx2
!endif 

indx=ix2+(nx-nx2)/2


!if (iy2< (ny2-1)/2) then
!indy=iy2
!else
!indy=iy2+ny-ny2
!endif 

indy=iy2+(ny-ny2)/2


if (iz2< (nz2-1)/2) then
indz=iz2
else
indz=iz2+nz-nz2
endif 

indxx2=indx + indy*nx + indz*nx*ny		
!print *, ix2,nx2,iy2,ny2,iz2,nz2,nx,ny,indxx2

end function indxx2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function indzz2(iz2,nz2)
integer, intent(in)	::iz2,nz2
integer	::indz


if (iz2<= (nz2-1)/2) then
indz=iz2
else
indz=iz2+nz-nz2
endif 

indzz2=indz

end function indzz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_kr(e_kr_tot,w_kr_tot)
real(kind(idp)), allocatable	:: SF_kr(:,:,:), array_c(:,:), array_f(:,:)
complex(KIND=idpc),intent(in) 	:: w_kr_tot(0:size_h_kr-1,0:size_h_kr-1,0:nk2d-1)
real (kind=idp),intent(in)	:: e_kr_tot(0:size_h_kr-1,0:nk2d-1)
 CHARACTER(2) :: alpha


!from SCF eigenvalues, z but not k resolved

allocate(SF_kr(0:nz-1, 0:1, 0:nk2d*size_h_kr-1))


SF_kr=0

do ik=0, nk2d-1			!kpoint
 do i=0, size_h_kr-1		!eigenvalue 
  do iz=0,nz-1			!layer
   do icf=0,1			!c/f
    do is=0,1			!up/down sum) 

     indp=iz + icf*nz + is*2*nz
     ind2p=ik+i*nk2d
     sf_kr(iz,icf,ind2p)=sf_kr(iz,icf,ind2p)+abs(w_kr_tot(indp,i,ik))**2*k_vecs2d(ik,4)/real(nkx*nky)
    enddo
   enddo
  enddo
 enddo
enddo
!!!!!!!!


do iz=0,nz-1
 WRITE(UNIT=alpha,FMT='(I2.2)') iz
 OPEN(unit=100+iz,file=trim(label)//'/sf/sf_kr_scf_nz'//alpha,status='unknown')



  write(100+iz,"(i10)") size_h_kr*nk2d

 do ik=0, nk2d-1			!kpoint

  allocate(array_c(0:size_h_kr-1,2))
  allocate(array_f(0:size_h_kr-1,2))
  do i=0, size_h_kr-1		!eigenvalue 

array_c(i,1)=e_kr_tot(i,ik)
array_c(i,2)=sf_kr(iz,0,ind2p)
array_f(i,1)=e_kr_tot(i,ik)
array_f(i,2)=sf_kr(iz,1,ind2p)

!call sort_and_squeeze(array_c)

     ind2p=ik+i*nk2d

     write(100+iz,format_smear) e_kr_tot(i,ik), sf_kr(iz,0,ind2p), sf_kr(iz,1,ind2p)
 !    write(100+iz,format_smear) array_c(i,1), array_c(i,2), array_f(i,2)

  enddo
  deallocate (array_c)
  deallocate (array_f)

 enddo
 

 close(unit=100+iz)
enddo

!not resolved (average over z)

 OPEN(unit=100,file=trim(label)//'/sf/sf_kr_scf',status='unknown')


  write(100,"(i10)") size_h_kr*nk2d

 do ik=0, nk2d-1			!kpoint
  do i=0, size_h_kr-1		!eigenvalue 

     ind2p=ik+i*nk2d

     write(100,format_smear) e_kr_tot(i,ik), sum(sf_kr(:,0,ind2p))/nz, sum(sf_kr(:,1,ind2p))/nz

  enddo
 enddo
 
 close(unit=100)



deallocate(sf_kr)


end subroutine spectralfunction_kr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_kr_nscf(e_kr_tot_nscf,w_kr_tot_nscf)
real(kind(idp)), allocatable	:: SF_kr_nscf(:,:,:,:),SF_kr_nscf_tot(:,:,:)
complex(KIND=idpc),intent(in) 	:: w_kr_tot_nscf(0:size_h_kr_red-1,0:size_h_kr_red_e-1,0:nkr_plot-1)
real (kind=idp),intent(in)	:: e_kr_tot_nscf(0:size_h_kr-1,0:nkr_plot-1)
 CHARACTER(2) :: alpha
 CHARACTER(3) :: beta
 CHARACTER(5) :: ak
 integer:: unit_file

!arpes, energy resolved, along path
if (lpathnscf) then

 allocate(SF_kr_nscf(0:nz-1, 0:1, 0:size_h_kr-1, 0:nkr_plot-1))
 SF_kr_nscf=0

 do ik=0, nkr_plot-1		!kpoint
  do i=0, size_h_kr-1		!eigenvalue 
   do iz=0,nz-1			!layer
    do icf=0,1			!c/f
     do is=0,1			!up/down sum) 

     indp=iz + icf*nz + is*2*nz
     sf_kr_nscf(iz,icf,i,ik)=sf_kr_nscf(iz,icf,i,ik)+abs(w_kr_tot_nscf(indp,i,ik))**2

     enddo
    enddo
   enddo
  enddo
 enddo


!! k and z resolved

do iz=0,nz-1
 do ik=0, nkr_plot-1
 
 unit_file=100+iz+ik*nz
 
  WRITE(UNIT=alpha,FMT='(I2.2)') iz
  WRITE(UNIT=beta,FMT='(I3.3)') ik
  
  OPEN(unit=unit_file,file=trim(label)//'/sf/sf_kr_nscf_nz'//alpha//'_k'//beta,status='unknown')


  write(unit_file,"(i10)") size_h_kr

  do i=0, size_h_kr-1		!eigenvalue 
 
   write(unit_file,format_smear) e_kr_tot_nscf(i,ik), sf_kr_nscf(iz,0,i,ik), sf_kr_nscf(iz,1,i,ik)
 
  enddo
 
  close(unit=unit_file)
 enddo
enddo

!k resolved, average over z
 do ik=0, nkr_plot-1
 
 unit_file=100+ik
 
  WRITE(UNIT=beta,FMT='(I3.3)') ik
  
  OPEN(unit=unit_file,file=trim(label)//'/sf/sf_kr_nscf_k'//beta,status='unknown')


  write(unit_file,"(i10)") size_h_kr

  do i=0, size_h_kr-1		!eigenvalue 
 
   write(unit_file,format_smear) e_kr_tot_nscf(i,ik), sum(sf_kr_nscf(:,0,i,ik))/real(nz), sum(sf_kr_nscf(:,1,i,ik))/real(nz)
 
  enddo
 
  close(unit=unit_file)
 enddo


deallocate(sf_kr_nscf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!on all k points - not lpathnscf


else

allocate(SF_kr_nscf_tot(0:nz_plot-1, 0:1, 0:nkr_plot*size_h_kr-1))


SF_kr_nscf_tot=0

do ik=0, nkr_plot-1		!kpoint
 do i=0, size_h_kr-1		!eigenvalue 
  do iz_ind=0,nz_plot-1	!z
   iz=z_values_plot(iz_ind)
   do icf=0,1			!c/f
    do is=0,1			!up/down sum) 

     indp=iz + icf*nz + is*2*nz
     ind2p=ik+i*nkr_plot
     sf_kr_nscf_tot(iz_ind,icf,ind2p)=sf_kr_nscf_tot(iz_ind,icf,ind2p)+abs(w_kr_tot_nscf(indp,i,ik))**2*kr_vecs_nscf(ik,4)/nkx_plot/nky_plot
     
    enddo
   enddo
  enddo
 enddo
enddo
!!!!!!!!
!sum over k (for dos)

do iz_ind=0,nz_plot-1	!z
 iz=z_values_plot(iz_ind)

 WRITE(UNIT=alpha,FMT='(I2.2)') iz
 unit_file=100
 OPEN(unit=unit_file,file=trim(label)//'/sf/sf_kr_nscf_nz'//alpha,status='unknown')

  write(unit_file,"(i10)") size_h_kr*nkr_plot

 do ik=0, nkr_plot-1			!kpoint
  do i=0, size_h_kr-1		!eigenvalue 

     ind2p=ik+i*nkr_plot

     write(unit_file,format_smear) e_kr_tot_nscf(i,ik), sf_kr_nscf_tot(iz_ind,0,ind2p), sf_kr_nscf_tot(iz_ind,1,ind2p)

  enddo
 enddo
 
 close(unit=unit_file)
enddo


!average over z and K (for total dos)
if (nz_plot==nz) then
 unit_file=100
 
 OPEN(unit=unit_file,file=trim(label)//'/sf/sf_kr_nscf',status='unknown')
  
  write(unit_file,"(i10)") size_h_kr*nkr_plot

 do ik=0, nkr_plot-1			!kpoint
  do i=0, size_h_kr-1		!eigenvalue 

     ind2p=ik+i*nkr_plot

     write(unit_file,format_smear) e_kr_tot_nscf(i,ik), sum(sf_kr_nscf_tot(:,0,ind2p))/nz, sum(sf_kr_nscf_tot(:,1,ind2p))/nz

  enddo
 enddo
 
 close(unit=unit_file)
endif
!!!!!!!!!!!!!!!!!!!!!!!!
!k,z resolved (arpes)
if(nz_plot<nz) then

do iz_ind=0,nz_plot-1	!z
 iz=z_values_plot(iz_ind)
 do ik=0, nkr_plot-1			!kpoint

 WRITE(UNIT=ak,FMT='(I5.5)') ik
 WRITE(UNIT=alpha,FMT='(I2.2)') iz
 unit_file=100
 
 OPEN(unit=unit_file,file=trim(label)//'/sf/sf_kr_arpes_nz'//alpha//'_k'//ak,status='unknown')

  
  write(unit_file,"(i10)") size_h_kr

  do i=0, size_h_kr-1		!eigenvalue 
     ind2p=ik+i*nkr_plot
    write(unit_file,format_smear) e_kr_tot_nscf(i,ik), sf_kr_nscf_tot(iz_ind,0,ind2p), sf_kr_nscf_tot(iz_ind,1,ind2p)
  enddo
    write(unit_file,format_smear) kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2),kr_vecs_nscf(ik,4)

 close(unit=unit_file)

 enddo
enddo

endif

deallocate(sf_kr_nscf_tot)


endif


end subroutine spectralfunction_kr_nscf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_r(e_r_tot,w_r_tot)
real (kind=idp)		:: e_r_tot(0:total_size-1,0:nk_r-1), occup
complex(KIND=idpc) 	:: w_r_tot(0:total_size-1,0:total_size-1,0:nk_r-1)
real(kind(ids)), allocatable	:: SF_r(:,:,:,:,:,:)
real(kind(ids)), allocatable	:: SF_r2(:,:,:,:,:,:)
real(kind(ids)), allocatable	:: SF_binning(:,:,:,:,:), E_binning(:)
real(kind(ids)), allocatable	:: SF_binning2(:,:,:,:,:)
real (kind=idp)		:: step, emin, emax

 CHARACTER(2) :: ax, ay, az
 CHARACTER(3) :: ak, ak2
integer	::unitfile, ik, ik2
real (kind=idp)	::weight_k
integer				::ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8
integer	::N_intervals, index_e, iz_ind
!integer, intent(in) ::nk_r

!from SCF eigenvalues, site resolved


print *, "Doing sf"

if (lscfr .or. .not. lpathnscfr) then
print *, "sf_r"

if (lbinning) then

 step=T/10
 emin=2*tc_
 emax=-emin


 n_intervals=int((emax-emin)/step)
 allocate(SF_binning(0:nx-1,0:ny-1,0:nz_plot-1, 0:2, 0:n_intervals))
 allocate(E_binning(0:n_intervals))

 SF_binning=0
 E_binning=0

 do j=0, n_intervals
  E_binning(j)=emin+j*step
 enddo


do ik=0,nk_r-1
 weight_k=k_vecs_r(ik,4)/nk_r_tot
 do i=0, total_size-1		!eigenvalue 
 
  if(e_r_tot(i,ik)>emin .and. e_r_tot(i,ik)<emax) then
   index_e=int((e_r_tot(i,ik)-emin)/step)

   do ix=0,nx-1	!x
    do iy=0,ny-1	!y
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
      do icf=0,1	!c/f
       do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind =ix + iy*nx + iz*nx*ny 

if (icf==0)     sf_binning(ix,iy,iz,icf,index_e)=sf_binning(ix,iy,iz,icf,index_e)+abs(w_r_tot(indp,i,ik))**2*weight_k
if (icf==1)     sf_binning(ix,iy,iz,icf,index_e)=sf_binning(ix,iy,iz,icf,index_e)+abs(w_r_tot(indp,i,ik))**2*weight_k!*b(ind)**2

icf2=1-icf

do is2=0,1

     ind2p=ix + iy*nx + iz*nx*ny + icf2*nx*ny*nz + is2*2*nx*ny*nz

     sf_binning(ix,iy,iz,2,index_e)=sf_binning(ix,iy,iz,2,index_e)+real(conjg(w_r_tot(indp,i,ik))*(w_r_tot(ind2p,i,ik)))*weight_k!*b(ind)

enddo
     

     
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 enddo
enddo



if (lkred_r .and. nx==ny) then

allocate(SF_binning2(0:nx-1,0:ny-1,0:nz_plot-1, 0:1, 0:n_intervals))


 do ix=0,nx-1	!x
  do iy=0,ny-1	!y

sf_binning2(ix,iy,:,:,:)=(sf_binning(ix,iy,:,:,:)+&
			 sf_binning(nx-1-ix,iy,:,:,:)+&
			 sf_binning(ix,ny-1-iy,:,:,:)+&
			 sf_binning(nx-1-ix,ny-1-iy,:,:,:)+&
			 sf_binning(iy,ix,:,:,:)+&
			 sf_binning(ny-1-iy,ix,:,:,:)+&
			 sf_binning(iy,nx-1-ix,:,:,:)+&
			 sf_binning(ny-1-iy,nx-1-ix,:,:,:))/8

  enddo
 enddo

sf_binning=sf_binning2

deallocate(sf_binning2)

endif



!!!!!
 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)

    WRITE(UNIT=ax,FMT='(I2.2)') ix
    WRITE(UNIT=ay,FMT='(I2.2)') iy
    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    ind=ix +iy*nx +iz*nx*ny
    unitfile=100+ind
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nx'//ax//'_ny'//ay//'_nz'//az,status='unknown')

    write(unitfile,"(i10)") n_intervals+1

   do j=0, n_intervals
     write(unitfile,format_smear) e_binning(j), sf_binning(ix,iy,iz,0,j), sf_binning(ix,iy,iz,1,j), sf_binning(ix,iy,iz,2,j)
   enddo

     write(unitfile,format_smear) b(ind)

 
    close(unit=unitfile)
 
   enddo
  enddo
 enddo






deallocate(SF_binning)



else




allocate(SF_r(0:nx-1,0:ny-1,0:nz_plot-1, 0:2, 0:total_size-1,0:nk_r-1))


!0:c
!1:f
!2:cf

SF_r=0


do ik=0,nk_r-1
 weight_k=k_vecs_r(ik,4)/nk_r_tot
 do i=0, total_size-1		!eigenvalue 
  do ix=0,nx-1	!x
   do iy=0,ny-1	!y
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
     do icf=0,1	!c/f
      do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind =ix + iy*nx + iz*nx*ny 

if (icf==0)     sf_r(ix,iy,iz,icf,i,ik)=sf_r(ix,iy,iz,icf,i,ik)+abs(w_r_tot(indp,i,ik))**2*weight_k
if (icf==1)     sf_r(ix,iy,iz,icf,i,ik)=sf_r(ix,iy,iz,icf,i,ik)+abs(w_r_tot(indp,i,ik))**2*weight_k*b(ind)**2	!hopping into the tip

icf2=1-icf

do is2=0,1

     ind2p=ix + iy*nx + iz*nx*ny + icf2*nx*ny*nz + is2*2*nx*ny*nz

     sf_r(ix,iy,iz,2,i,ik)=sf_r(ix,iy,iz,2,i,ik)+real(conjg(w_r_tot(indp,i,ik))*(w_r_tot(ind2p,i,ik)))*weight_k*b(ind)

enddo
     

     
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
!!!!!!!!

!correction for reduced k set

!if (lkred_r .and. nx==ny) then

!allocate(SF_r2(0:nx-1,0:ny-1,0:nz-1, 0:1, 0:total_size-1,0:nk_r-1))


! do ix=0,nx-1	!x
!  do iy=0,ny-1	!y

!sf_r2(ix,iy,:,:,:,:)=	(sf_r(ix,iy,:,:,:,:)+&
!			 sf_r(nx-1-ix,iy,:,:,:,:)+&
!			 sf_r(ix,ny-1-iy,:,:,:,:)+&
!			 sf_r(nx-1-ix,ny-1-iy,:,:,:,:)+&
!			 sf_r(iy,ix,:,:,:,:)+&
!			 sf_r(ny-1-iy,ix,:,:,:,:)+&
!			 sf_r(iy,nx-1-ix,:,:,:,:)+&
!			 sf_r(ny-1-iy,nx-1-ix,:,:,:,:))/8


!  enddo
! enddo

!sf_r=sf_r2

!deallocate(sf_r2)

!endif


!!!!!
!write
 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

    WRITE(UNIT=ax,FMT='(I2.2)') ix
    WRITE(UNIT=ay,FMT='(I2.2)') iy
    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    unitfile=100+ix +iy*nx +iz*nx*ny
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nx'//ax//'_ny'//ay//'_nz'//az,status='unknown')

    write(unitfile,"(i10)") total_size*nk_r

   do ik=0,nk_r-1
    do i=0, total_size-1		!eigenvalue 
     write(unitfile,format_smear) e_r_tot(i,ik), sf_r(ix,iy,iz,0,i,ik), sf_r(ix,iy,iz,1,i,ik), sf_r(ix,iy,iz,2,i,ik)
    enddo
   enddo
 
    close(unit=unitfile)
 
   enddo
  enddo
 enddo

!average over x and y

     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    unitfile=100+iz
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nz'//az,status='unknown')

    write(unitfile,"(i10)") total_size*nk_r

   do ik=0,nk_r-1
    do i=0, total_size-1		!eigenvalue 
     write(unitfile,format_smear) e_r_tot(i,ik), sum(sf_r(:,:,iz,0,i,ik))/(real(nx*ny)), sum(sf_r(:,:,iz,1,i,ik))/(real(nx*ny)),sum(sf_r(:,:,iz,2,i,ik))/(real(nx*ny))
    enddo
   enddo
    
    close(unit=unitfile)
 
   enddo



!average over x, y, z


    unitfile=100
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r',status='unknown')

    write(unitfile,"(i10)") total_size*nk_r

   do ik=0,nk_r-1
    do i=0, total_size-1		!eigenvalue 
     write(unitfile,format_smear) e_r_tot(i,ik), sum(sf_r(:,:,:,0,i,ik))/(real(nx*ny*nz)), sum(sf_r(:,:,:,1,i,ik))/(real(nz*ny*nz)),sum(sf_r(:,:,:,2,i,ik))/(real(nz*ny*nz))
    enddo
   enddo
   
    close(unit=unitfile)

deallocate(sf_r)

endif

endif

!!!!!!!!!!!as a function of external momentum, for lpathnscfr

!nkpath_r=nkpath_r+1
!nk_r=3*nkpath_r
!nk_r_rot=(nx-1)/2
!n_k_total=3*(nkpath_r/2+nkpath_r*nk_r_rot+1)

!do ind=0, n_k_total-1
!k_tot(ind,1)=(k_vecs_r(list_kk(ind,3),1)+list_kk(ind,1))/nx
!k_tot(ind,2)=(k_vecs_r(list_kk(ind,3),2)+list_kk(ind,2))/nx
!write(*, "(i5, 2f10.5)") ind,  k_tot(ind,1), k_tot(ind,2)
!enddo


print *, "End sf"

end subroutine spectralfunction_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kr_path(kr_vecs, nkpath_plot, nkr_plot)
real(kind=idp),intent(out):: kr_vecs(0:nkr_plot-1,4)
integer, intent(in)	::nkpath_plot, nkr_plot
integer	::ikx

!G-X-M-G
if (nx_kr==1 .and. ny_kr==1 .and. (surf_type=="100" .or. surf_type=="001") .and. .true.) then

  do ikx=0,nkxpath_plot-1
  
   kr_vecs (ikx,1)=real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx,2)=0

  enddo

  do ikx=0,nkypath_plot-1

   kr_vecs (nkxpath_plot+ikx,1)=0.5
   kr_vecs (nkxpath_plot+ikx,2)=real(ikx)/real(2*nkypath_plot)

  enddo

 do ikx=0,nkpath_plot_xy-1
 
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,1)=0.5-real(ikx)/real(2*nkpath_plot_xy)
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,2)=kr_vecs (nkxpath_plot+nkypath_plot+ikx,1)

!   kr_vecs (3*nkpath_plot+nkpath_plot_xy+ikx,1)=-real(ikx)/real(2*nkpath_plot_xy)
!   kr_vecs (3*nkpath_plot+nkpath_plot_xy+ikx,2)=0.5-real(ikx)/real(2*nkpath_plot_xy)

 enddo


!M-G-X-M
elseif (nx_kr==1 .and. ny_kr==1 .and. (surf_type=="100" .or. surf_type=="001") .and. .false.) then


 do ikx=0,nkpath_plot_xy-1
 
   kr_vecs (ikx,1)=0.5-real(ikx)/real(2*nkpath_plot_xy)
   kr_vecs (ikx,2)=kr_vecs (ikx,1)

 enddo

  do ikx=0,nkxpath_plot-1
  
   kr_vecs (nkpath_plot_xy+ikx,1)=real(ikx)/real(2*nkxpath_plot)
   kr_vecs (nkpath_plot_xy+ikx,2)=0

  enddo

  do ikx=0,nkypath_plot

   kr_vecs (nkxpath_plot+nkpath_plot_xy+ikx,1)=0.5
   kr_vecs (nkxpath_plot+nkpath_plot_xy+ikx,2)=real(ikx)/real(2*nkypath_plot)

  enddo




!G-X-M-X'G but the BZ is here reduced
elseif (nx_kr==2 .and. ny_kr==1 .and. lfullbz  .and. (surf_type=="100" .or. surf_type=="001")) then

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx,1)=nx_kr*real(ikx)/real(4*nkxpath_plot)
   kr_vecs (ikx,2)=0
  enddo

  do ikx=0,nkypath_plot-1
   kr_vecs (nkxpath_plot+ikx,1)=0.5
   kr_vecs (nkxpath_plot+ikx,2)=real(ikx)/real(2*nkypath_plot)
  enddo

  do ikx=0,nkxpath_plot-1
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,1)=(0.5-real(ikx)/real(2*nkxpath_plot))
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,2)=0.5
  enddo

  do ikx=0,nkypath_plot-1
   kr_vecs (2*nkxpath_plot+nkypath_plot+ikx,1)=0
   kr_vecs (2*nkxpath_plot+nkypath_plot+ikx,2)=0.5-real(ikx)/real(2*nkypath_plot)
  enddo
 
elseif (nx_kr==2 .and. ny_kr==1 .and. .not. lfullbz  .and. (surf_type=="100" .or. surf_type=="001")) then

  do ikx=0,nkxpath_plot-1
  
   kr_vecs (ikx,1)=real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx,2)=0

  enddo

  do ikx=0,nkypath_plot-1

   kr_vecs (nkxpath_plot+ikx,1)=0.5
   kr_vecs (nkxpath_plot+ikx,2)=real(ikx)/real(2*nkypath_plot)

  enddo

  do ikx=0,nkxpath_plot-1
  
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,1)=0.5-real(ikx)/real(2*nkxpath_plot)
   kr_vecs (nkxpath_plot+nkypath_plot+ikx,2)=0.5

  enddo

  do ikx=0,nkypath_plot-1

   kr_vecs (2*nkxpath_plot+nkypath_plot+ikx,1)=0
   kr_vecs (2*nkxpath_plot+nkypath_plot+ikx,2)=0.5-real(ikx)/real(2*nkypath_plot)

  enddo



elseif(surf_type=="110") then
!GXSYG

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx,1)=real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx,2)=0
  enddo

  do ikx=0,nkypath_plot-1
   kr_vecs (ikx+nkxpath_plot,1)=0.5
   kr_vecs (ikx+nkxpath_plot,2)=real(ikx)/real(2*nkypath_plot)
  enddo

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx+nkxpath_plot+nkypath_plot,1)=0.5-real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx+nkxpath_plot+nkypath_plot,2)=0.5
  enddo

  do ikx=0,nkypath_plot
   kr_vecs (ikx+2*nkxpath_plot+nkypath_plot,1)=0
   kr_vecs (ikx+2*nkxpath_plot+nkypath_plot,2)=0.5-real(ikx)/real(2*nkypath_plot)
  enddo


elseif(surf_type=="210") then
!SYGXS so that gamma is at the center

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx,1)=0.5-real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx,2)=0.5
  enddo

  do ikx=0,nkypath_plot-1
   kr_vecs (ikx+nkxpath_plot,1)=0
   kr_vecs (ikx+nkxpath_plot,2)=0.5-real(ikx)/real(2*nkypath_plot)
  enddo

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx+nkxpath_plot+nkypath_plot,1)=real(ikx)/real(2*nkxpath_plot)
   kr_vecs (ikx+nkxpath_plot+nkypath_plot,2)=0
  enddo

  do ikx=0,nkypath_plot
   kr_vecs (ikx+2*nkxpath_plot+nkypath_plot,1)=0.5
   kr_vecs (ikx+2*nkxpath_plot+nkypath_plot,2)=real(ikx)/real(2*nkypath_plot)
  enddo



elseif(surf_type=="111") then
!GMKG
!G=0,0	M=1/2,1/2	K=2/3,1/3


  do ikx=0,nkypath_plot-1
   kr_vecs (ikx,1)=real(ikx)/real(2*nkypath_plot)
   kr_vecs (ikx,2)=kr_vecs(ikx,1)
  enddo

  do ikx=0,nkxpath_plot-1
   kr_vecs (ikx+nkypath_plot,1)=0.5+real(ikx)/real(nkxpath_plot)/6.0d0
   kr_vecs (ikx+nkypath_plot,2)=0.5-real(ikx)/real(nkxpath_plot)/6.0d0
  enddo

  do ikx=0,nkxypath_plot
   kr_vecs (ikx+nkypath_plot+nkxpath_plot,1)=2./3.*(1-real(ikx)/real(nkxypath_plot))
   kr_vecs (ikx+nkypath_plot+nkxpath_plot,2)=1./3.*(1-real(ikx)/real(nkxypath_plot))
  enddo



endif


! kr_vecs_nscf (3*nkpath_plot+ikx,1)=real(ikx)/real(2*nkpath_plot)
! kr_vecs_nscf (3*nkpath_plot+ikx,2)=-real(ikx)/real(2*nkpath_plot)

! kr_vecs_nscf (4*nkpath_plot+ikx,1)=0.5-real(ikx)/real(2*nkpath_plot)
! kr_vecs_nscf (4*nkpath_plot+ikx,2)=-0.5

! kr_vecs_nscf (5*nkpath_plot+ikx,1)=0
! kr_vecs_nscf (5*nkpath_plot+ikx,2)=-0.5+real(ikx)/real(2*nkpath_plot)

  
!   kr_vecs (4*nkpath_plot+2*nkpath_plot_xy,1)=0
!   kr_vecs (4*nkpath_plot+2*nkpath_plot_xy,2)=0
!   kr_vecs (2*nkpath_plot+nkpath_plot_xy,1)=0
!   kr_vecs (2*nkpath_plot+nkpath_plot_xy,2)=0



kr_vecs(:,3)=0		!z
kr_vecs(:,4)=1		!weight



end subroutine set_kr_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kr_path_r(k_vecs_r,nkpath_r,nk_r)
real(kind=idp),intent(out), allocatable:: k_vecs_r(:,:)
integer,intent(inout)	::nkpath_r
integer,intent(out)	::nk_r
integer	::ikx

!reduced procedure

if (nx==ny .and. lsfrkred) then

nkpath_x=nkpath_r
nkpath_y=nkpath_x

!print *, nkpath_r, nkpath_r_xy


!!!!!!!!!!!!!!!!
!k in the reduced zone

!nx=ny=odd
if (mod(nx,2)==1) then

 nk_r=2*nkpath_r+nkpath_r_xy+1
 allocate(k_vecs_r(0:nk_r-1,4))


  do ikx=0,nkpath_r-1
   k_vecs_r (ikx,1)=real(ikx)/real(2*nkpath_r)
   k_vecs_r (ikx,2)=0
  enddo
  

  do ikx=0,nkpath_r-1
   k_vecs_r (ikx+nkpath_r,1)=0.5
   k_vecs_r (ikx+nkpath_r,2)=real(ikx)/real(2*nkpath_r)
  enddo
  
  
  do ikx=0,nkpath_r_xy-1
   k_vecs_r (ikx+nkpath_r*2,1)=0.5-real(ikx)/real(2*nkpath_r_xy)
   k_vecs_r (ikx+nkpath_r*2,2)=0.5-real(ikx)/real(2*nkpath_r_xy)
  enddo

!nx=ny=even
 
else

 nk_r=2*(nkpath_r+1)+nkpath_r_xy+1
 allocate(k_vecs_r(0:nk_r-1,4))

  
  do ikx=0,nkpath_r
   k_vecs_r (ikx,1)=real(ikx)/real(2*nkpath_r)
   k_vecs_r (ikx,2)=0
  enddo
  

  do ikx=0,nkpath_r
   k_vecs_r (ikx+nkpath_r+1,1)=0
   k_vecs_r (ikx+nkpath_r+1,2)=real(ikx)/real(2*nkpath_r)
  enddo
  
  
  do ikx=0,nkpath_r_xy
   k_vecs_r (ikx+nkpath_r*2+2,1)=0.5-real(ikx)/real(2*nkpath_r_xy)
   k_vecs_r (ikx+nkpath_r*2+2,2)=0.5-real(ikx)/real(2*nkpath_r_xy)
  enddo
  
endif

!gamma
k_vecs_r(nk_r-1,1)=0
k_vecs_r(nk_r-1,2)=0


k_vecs_r(:,3)=0		!z
k_vecs_r(:,4)=1		!weight



!!!!!!!!!!!!!!!!!!!!!!!
!combined k-K (in the full BZ)

if (mod(nx,2)==1) nk_r_rot=(nx-1)/2		!5: -2 -1 0 +1 +2 --> 2
if (mod(nx,2)==0) nk_r_rot=nx/2			!4: -2 -1 0 +1 +2 --> 2

!if (mod(nx,2)==1) n_k_total=(2*nkpath_r+nkpath_r_xy)*(2*nk_r_rot+1)+1
!if (mod(nx,2)==0) n_k_total=(2*nkpath_r+nkpath_r_xy)*(2*nk_r_rot)+1

n_k_total=(2*nkpath_r+nkpath_r_xy)*nx+1

print *, "nk_r_rot, nkpath_r, nkpath_r_xy, n_k_total, nk_r"
print *, nk_r_rot, nkpath_r,nkpath_r_xy,n_k_total, nk_r

if (allocated(list_kk)) deallocate(list_kk)
allocate(list_kk(0:n_k_total-1,5))
allocate(k_tot(0:n_k_total-1,2))



!!!!!!!nx odd
if (mod(nx,2)==1) then

!along x

ind=-1

 do j=0, nkpath_r	! external
  ind=ind+1
 list_kk(ind,1)=0		!internal momentum x	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=j	!index of external momentum	(small)
 list_kk(ind,4)=1	! plus or minus external momentum
 list_kk(ind,5)=1	! plus or minus external momentum
!print *, ind,list_kk(ind,3)
 enddo


 do i=1, nk_r_rot

  do j=1, nkpath_r	! external
   ind=ind+1
 list_kk(ind,1)=i		!internal momentum	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=nkpath_r-j	!external momentum	(small)
 list_kk(ind,4)=-1	! plus or minus external momentum
 list_kk(ind,5)=1	!plus or minus external momentum
!print *, ind,list_kk(ind,3)
  enddo
  do j=1, nkpath_r	! external
   ind=ind+1
 list_kk(ind,1)=i		!internal momentum	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=j	!external momentum	(small)
 list_kk(ind,4)=1	! plus or minus external momentum
 list_kk(ind,5)=1	
!print *, "x", ind,list_kk(ind,3)
 enddo

enddo


!along y

 do j=0, nkpath_r-1	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=0		!internal momentum	(large)
list_kk(ind,3)=j+nkpath_r+1	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *,"xy", ind,list_kk(ind,3)
enddo

do i=1, nk_r_rot

 do j=1, nkpath_r	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=2*nkpath_r-j	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo
 do j=1, nkpath_r	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=nkpath_r+j	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
 enddo

enddo

!along x-y
do i=0, nk_r_rot-1
 do j=0, nkpath_r_xy-1	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,2)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,3)=j+2*nkpath_r+1	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
enddo
 do j=0, nkpath_r_xy-1	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,2)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,3)=nkpath_r_xy-j+2*nkpath_r-1	!external momentum	(small)
list_kk(ind,4)=-1	! plus or minus external momentum
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo

enddo

 do j=0, nkpath_r_xy-1	! external
  ind=ind+1
list_kk(ind,1)=0		!internal momentum	(large)
list_kk(ind,2)=0		!internal momentum	(large)
list_kk(ind,3)=j+2*nkpath_r+1	!external momentum	(small)
list_kk(ind,4)=1	
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
enddo

!!!!!!!!!!!!!
!nx=ny even

else if (mod(nx,2)==0) then

ind=-1

i=0
 do j=0, nkpath_r	! external
  ind=ind+1
 list_kk(ind,1)=i		!internal momentum x	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=j	!index of external momentum	(small)
 list_kk(ind,4)=1	! plus or minus external momentum
 list_kk(ind,5)=1	! plus or minus external momentum
!print *, ind,list_kk(ind,3)
 enddo


 do i=1, nk_r_rot-1

  do j=1, nkpath_r	! external
   ind=ind+1
 list_kk(ind,1)=i		!internal momentum	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=nkpath_r-j	!external momentum	(small)
 list_kk(ind,4)=-1	! plus or minus external momentum
 list_kk(ind,5)=1	!plus or minus external momentum
!print *, ind,list_kk(ind,3)
  enddo
  do j=1, nkpath_r	! external
   ind=ind+1
 list_kk(ind,1)=i		!internal momentum	(large)
 list_kk(ind,2)=0		!internal momentum y	(large)
 list_kk(ind,3)=j	!external momentum	(small)
 list_kk(ind,4)=1	! plus or minus external momentum
 list_kk(ind,5)=1	
!print *, "x", ind,list_kk(ind,3)
 enddo

enddo

i=nk_r_rot
 do j=1, nkpath_r-1	! external
  ind=ind+1
list_kk(ind,1)=i		!internal momentum x	(large)
list_kk(ind,2)=0		!internal momentum y	(large)
list_kk(ind,3)=nkpath_r-j	!index of external momentum	(small)
list_kk(ind,4)=-1	
list_kk(ind,5)=1	
!print *, "x", ind,list_kk(ind,3)
enddo


!along y

i=0
 do j=0, nkpath_r	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=j+nkpath_r+1	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *,"y", ind,list_kk(ind,3)
enddo

do i=1, nk_r_rot-1

 do j=1, nkpath_r	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=2*nkpath_r+1-j	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo
 do j=1, nkpath_r	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=nkpath_r+1+j	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
 enddo

enddo

  i=nk_r_rot
 do j=1, nkpath_r-1	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot		!internal momentum	(large)
list_kk(ind,2)=i		!internal momentum	(large)
list_kk(ind,3)=2*nkpath_r+1-j	!external momentum	(small)
list_kk(ind,4)=1	
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo


!along x-y

  i=0
do j=0,nkpath_r_xy	! external
 ind=ind+1
list_kk(ind,1)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,2)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,3)=nkpath_r_xy-j+2*nkpath_r+2	!external momentum	(small)
list_kk(ind,4)=-1	
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo


do i=1, nk_r_rot-1
 do j=1, nkpath_r_xy	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,2)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,3)=j+2*nkpath_r+2	!external momentum	(small)
list_kk(ind,4)=1	! plus or minus external momentum
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
enddo
 do j=1, nkpath_r_xy	! external
  ind=ind+1
list_kk(ind,1)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,2)=nk_r_rot-i		!internal momentum	(large)
list_kk(ind,3)=nkpath_r_xy-j+2*nkpath_r+2	!external momentum	(small)
list_kk(ind,4)=-1	! plus or minus external momentum
list_kk(ind,5)=-1	
!print *, ind,list_kk(ind,3)
enddo

enddo

 do j=1, nkpath_r_xy	! external
  ind=ind+1
list_kk(ind,1)=0		!internal momentum	(large)
list_kk(ind,2)=0		!internal momentum	(large)
list_kk(ind,3)=j+2*nkpath_r+2	!external momentum	(small)
list_kk(ind,4)=1	
list_kk(ind,5)=1	
!print *, ind,list_kk(ind,3)
enddo


endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *,ind,n_k_total

do ik=0, n_k_total-1
!print *, ik
k_tot(ik,1)=(k_vecs_r(list_kk(ik,3),1)+list_kk(ik,4)*list_kk(ik,1))/nx
k_tot(ik,2)=(k_vecs_r(list_kk(ik,3),2)+list_kk(ik,5)*list_kk(ik,2))/nx
!write(*, "(2i5, 2f10.5)") ik,  k_tot(ik,1), k_tot(ik,2)
enddo


!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else !nx .ne. ny, or lsfrkred=false
!k in the full bz, or in the reducd one, no k+K procedure
!stop

!if (mod(nkpath_plot,2)==1) nkpath_plot=nkpath_plot-1	!it must be even to contain gamma (in the reduced zone)

print *, "lfullbz", lfullbz

if (lfullbz) then 
nxn=nx
nyn=ny
else
nxn=1
nyn=1
endif

!print *, nxn, nyn

!nkpath_plot_xy=int(sqrt(2.0D0)*nkpath_plot)


 nk_r=2*nkpath_plot+nkpath_plot_xy+1
 n_k_total=nk_r

if(allocated(list_kk)) deallocate(list_kk)
allocate(k_vecs_r(0:nk_r-1,4))
allocate(list_kk(0:n_k_total-1,4))
allocate(k_tot(0:n_k_total-1,2))


!!!!!!!!!!!!!!!!
!k in the reduced (or full) zone


!along x
  do ikx=0,nkpath_plot-1
   k_vecs_r (ikx,1)=real(ikx)/real(2*nkpath_plot)*nxn
   k_vecs_r (ikx,2)=0
  enddo
  

!  do ikx=0,2*nkpath_plot-1
!   k_vecs_r (ikx,1)=(-0.5+real(ikx)/real(2*nkpath_plot))*nxn
!   k_vecs_r (ikx,2)=0.25
!  enddo


!along y

  do ikx=0,nkpath_plot-1
   k_vecs_r (ikx+nkpath_plot,1)=0.5*nxn
   k_vecs_r (ikx+nkpath_plot,2)=real(ikx)/real(2*nkpath_plot)*nyn
  enddo
  
  !along xy
  
  do ikx=0,nkpath_plot_xy-1
   k_vecs_r (ikx+nkpath_plot*2,1)=(0.5-real(ikx)/real(2*nkpath_plot_xy))*nxn
   k_vecs_r (ikx+nkpath_plot*2,2)=(0.5-real(ikx)/real(2*nkpath_plot_xy))*nyn
  enddo


!  do ikx=0,nkpath_plot-1
!   k_vecs_r (ikx+nkpath_plot*3,1)=0
!   k_vecs_r (ikx+nkpath_plot*3,2)=(real(ikx)/real(2*nkpath_plot))*nyn
! enddo

 ! do ikx=0,nkpath_plot-1
 !  k_vecs_r (ikx+nkpath_plot*4,1)=(real(ikx)/real(2*nkpath_plot))*nxn
 !  k_vecs_r (ikx+nkpath_plot*4,2)=0.5*nyn
 ! enddo


!   k_vecs_r (nkpath_plot*5,1)=0.5*nxn
 !  k_vecs_r (nkpath_plot*5,2)=0.5*nyn

!gamma st the end

   k_vecs_r (nk_r-1,1)=0
   k_vecs_r (nk_r-1,2)=0



k_vecs_r(:,3)=0		!z
k_vecs_r(:,4)=1		!weight


!!!!!!!!!!!!!!!!!!!!!!!
!combined k-K

do ik=0, n_k_total-1
list_kk(ik,1)=0		!internal momentum x	(large)
list_kk(ik,2)=0		!internal momentum y	(large)
list_kk(ik,3)=ik	!index of external momentum	(small)
k_tot(ik,1)=k_vecs_r(ik,1)/nxn
k_tot(ik,2)=k_vecs_r(ik,2)/nyn
!write(*, "(3i5, 2f10.5, i5, 4f10.5)")   ik, k_tot(ik,1), k_tot(ik,2) 
enddo

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do ik=0, n_k_total-1

write(*, "(3i5, 2f10.5, i5, 4f10.5)")   ik, list_kk(ik,1),list_kk(ik,2),real(list_kk(ik,1))/nx,real(list_kk(ik,2))/ny,&
				list_kk(ik,3),k_vecs_r(list_kk(ik,3),1)/nx, k_vecs_r(list_kk(ik,3),2)/nx,&
				 k_tot(ik,1), k_tot(ik,2) 

enddo


!stop


!sottinteso che kx sia in unita 1/nx e ky in 1/ny


end subroutine set_kr_path_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_wrtot(e_r_tot, w_r_tot)
complex(KIND=idpc) ,intent(in)	:: w_r_tot(0:total_size-1,0:total_size-1, 0:nk_r-1)
real(KIND=idp) ,intent(in) 	:: E_r_tot(0:total_size-1, 0:nk_r-1)
 character(3)	::alpha
integer	::unit_file, unit_file_2


do ik=0, nk_r-1
 unit_file=100+ik
 WRITE(UNIT=alpha,FMT='(I3.3)') ik

 OPEN(unit=unit_file,file=trim(label)//"/wfc/wfc_r2_k"//alpha,status='unknown')

 do i=0, total_size-1			!eigenvalue
  do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z
     do icf=0,1
      do is=0,1    
     
      indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
      write (unit_file, format_wr) i, indp, real(w_r_tot(indp,i,ik)), aimag(w_r_tot(indp,i,ik))
      
      enddo     					
     enddo
    enddo
   enddo
  enddo
 enddo
 
 close(unit_file)

!eigenvalues
 unit_file_2=200+ik
 OPEN(unit=unit_file_2,file=trim(label)//"/wfc/e_r2_k"//alpha,status='unknown')

 do i=0, total_size-1			!eigenvalue
     
  write (unit_file_2, format_er) i, e_r_tot(i,ik)
     					
 enddo

 close(unit_file_2)

enddo

end subroutine write_wrtot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readwr(ik,e,w)
complex(KIND=idpc) 	:: W(0:total_size-1,0:total_size-1)
real(KIND=idp)		:: E(0:total_size-1), a, b
integer,intent(in)	:: ik
integer	::j
 character(3)	::alpha


print *, "Reading wfc..."


 WRITE(UNIT=alpha,FMT='(I3.3)') ik

 OPEN(unit=38,file=trim(label)//"/wfc/wfc_r2_k"//alpha,status='unknown')

 do i=0, total_size-1			!eigenvalue
  do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z
     do icf=0,1
      do is=0,1    
     
      indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
      read (38, format_wr) j, ind2p, a, b
      w(ind2p,j)=a+(0,1)*b     
      
      enddo     					
     enddo
    enddo
   enddo
  enddo
 enddo
 
 close(38)

!eigenvalues
 OPEN(unit=39,file=trim(label)//"/wfc/e_r2_k"//alpha,status='unknown')

 do i=0, total_size-1			!eigenvalue
     
     read (39, format_er) j, e(j)
     					
 enddo

 close(39)


print *, "... OK reading wfc"


end subroutine readwr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_k(e_k_tot,w_k_tot)
real(kind(idp)), allocatable	:: SF_k(:,:)
real (kind=idp)		:: e_k_tot(0:3,0:nk3d-1)
complex(KIND=idpc) 	:: w_k_tot(0:3,0:3,0:nk3d-1)


allocate(SF_k(0:1, 0:nk3d*4-1))


SF_k=0

do ik=0, nk3d-1			!kpoint
 do i=0, 3		!eigenvalue 
  do icf=0,1			!c/f
   do is=0,1			!up/down sum) 

     indp= icf*2 + is
     ind2p=ik+i*nk3d
     sf_k(icf,ind2p)=sf_k(icf,ind2p)+abs(w_k_tot(indp,i,ik))**2*k_vecs3d(ik,4)/real(nkx_k*nky_k*nkz)
     
   enddo
  enddo
 enddo
enddo
!!!!!!!!


 OPEN(unit=100,file=trim(label)//'/sf/sf_k_scf',status='unknown')


  write(100,"(i10)") 4*nk3d

 do ik=0, nk3d-1			!kpoint
  do i=0, 3		!eigenvalue 

     ind2p=ik+i*nk3d

     write(100,format_smear) e_k_tot(i,ik), sf_k(0,ind2p), sf_k(1,ind2p)

  enddo
 enddo
 
 close(unit=100)

deallocate(sf_k)


end subroutine spectralfunction_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_bands_r(e_r_tot, lscfr)
real (kind=idp)		:: e_r_tot(0:total_size-1,0:nk_r-1)
integer		::i,j,ik
logical, intent(in)	::lscfr

!print *, lscfr, lpathnscfr, lfullbz


if (lscfr) OPEN(unit=100,file=trim(label)//'/bands_r_scf',status='unknown')
if (.not. lscfr .and. .not. lpathnscfr) OPEN(unit=100,file=trim(label)//'/bands_r_nscf',status='unknown')
if (.not. lscfr .and. lpathnscfr .and. lfullbz) OPEN(unit=100,file=trim(label)//'/bands_r_path_full',status='unknown')
if (.not. lscfr .and. lpathnscfr .and. .not. lfullbz) OPEN(unit=100,file=trim(label)//'/bands_r_path_red',status='unknown')


do i=0,total_size-1
 do ik= 0, nk_r-1
  write(100,"(i5, 3f10.5, 5f15.10)") ik, k_vecs_r(ik,1),k_vecs_r(ik,2), k_vecs_r(ik,3),e_r_tot(i,ik)
 enddo
 write(100,"(a1)") ""
enddo
     
     close(unit=100)
     

end subroutine write_bands_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_index_site(ind, ix,iy,iz)
integer	::ind, ix,iy,iz

ix=mod(ind,nx)
iy=mod((ind-ix)/nx,ny)
iz=mod((ind-ix-iy*nx)/nx/ny,nz)

!ind=ix + iy*nx + iz*nx*ny



end subroutine convert_index_site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gcd(v, t)
  integer :: gcd
  integer, intent(in) :: v, t
  integer :: c, b, a
 
  b = t
  a = v
  do
     c = mod(a, b)
     if ( c == 0) exit
     a = b
     b = c
  end do
  gcd = b ! abs(b)
end function gcd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_r_k(ik)
integer	,intent(in)	::ik
real(kind(ids)), allocatable	:: SF_k(:,:,:)
complex(kind=ids),allocatable:: rot_kx(:,:,:), w_k(:,:)
 CHARACTER(2) :: ax, ay, az
 CHARACTER(3) :: ak, ak2
integer	::unitfile

!if (.not. lsfrkred .and. lfullbz) print *, "doing sfk"

if (ik==0) print *, n_k_total

if (ik==0) print *, "sfk"

!allocate(w_k(0:total_size-1,0:total_size-1))		!check size?
allocate(w_k(0:4*nz_plot-1,0:total_size-1))
allocate(sf_k(0:nz_plot-1,0:1,0:total_size-1))

if (ik==0 .and. lfullbz)  write(*, "(3a5, 2a10, a5,4a10)") "iktot", "nK_x","nK_y","K_x","K_y", "ik", "k_x", "k_y", "k_tot x", "k_tot y"

!do ik=0, n_k_total-1

if (lfullbz) write(*, "(3i5, 2f10.5, i5, 4f10.5)")   ik, list_kk(ik,1),list_kk(ik,2),real(list_kk(ik,1))/nx,real(list_kk(ik,2))/ny,&
				list_kk(ik,3),k_vecs_r(list_kk(ik,3),1)/nx, k_vecs_r(list_kk(ik,3),2)/nx,&
				 k_tot(ik,1), k_tot(ik,2) 


 
!eigenvectors in kx,ky,z basis
!print *, " rotation"

!print *, 4*nz, total_size

w_k=0

 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
!print *,  ix,iy, exp(2*(0,1)*pi*(list_kk(ik,1)*ix/real(nx)+list_kk(ik,2)*iy/real(ny)))
   do i=0,total_size-1
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
     do icf=0,1	!c/f
      do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind2p=iz + icf*nz + is*2*nz

if (.not. lsfrkred)  w_k(ind2p,i)= w_k(ind2p,i)+w(indp,i)                    * exp(-2*(0,1)*pi*( k_tot(ik,1)*ix/real(1)+k_tot(ik,2)*iy/real(1) ) )/sqrt(real(nx*ny))
if (lsfrkred)        w_k(ind2p,i)= w_k(ind2p,i)+w_r_tot(indp,i,list_kk(ik,3))* exp(-2*(0,1)*pi*( k_tot(ik,1)*ix/real(1)+k_tot(ik,2)*iy/real(1) ) )/sqrt(real(nx*ny))


       enddo 
      enddo
     enddo
    enddo
   enddo
  enddo


!spectral weights for k points
!print *, " sf_k"
sf_k=0

! do ix=0,nx-1	!x
!  do iy=0,ny-1	!y
   do i=0,total_size-1
!   do j=0,total_size-1
    do icf=0,1	!c/f
     do is=0,1	!up/down 
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

     ind2p=iz + icf*nz + is*2*nz
!     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz

      sf_k(iz,icf,i)=sf_k(iz,icf,i)+real((conjg(w_k(ind2p,i))*w_k(ind2p,i)))		!sum over spin
!     sf_k(iz,icf,i)=sf_k(iz,icf,i)+abs(conjg(w_k(indp,i))*w_k(indp,i))**2
!     sf_k(iz,icf,i)=sf_k(iz,icf,i)+abs(w_k0(ind2p,i))**2

       enddo
      enddo
     enddo
    enddo
!enddo
!enddo
!enddo
!print *, " print"

!sf_k=abs(sf_k)**2

!k resolved, average over z

    WRITE(UNIT=ak,FMT='(I3.3)') ik


    unitfile=100
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_k'//ak,status='unknown')

    write(unitfile,"(i10)") total_size

    do i=0, total_size-1		!eigenvalue 
if (lsfrkred)          write(unitfile,format_smear) e_r_tot(i,list_kk(ik,3)), sum(sf_k(:,0,i))/real(nz), sum(sf_k(:,1,i))/real(nz)
if (.not. lsfrkred)    write(unitfile,format_smear) e(i), sum(sf_k(:,0,i))/real(nz), sum(sf_k(:,1,i))/real(nz)
    enddo
    
   close(unit=unitfile)

!k-z resolved

     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    unitfile=100+iz
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nz'//az//'_k'//ak,status='unknown')

    write(unitfile,"(i10)") total_size
  
    do i=0, total_size-1		!eigenvalue 
    if (lsfrkred)  write(unitfile,format_smear) e_r_tot(i,list_kk(ik,3)), (sf_k(iz,0,i)), (sf_k(iz,1,i))
    if (.not. lsfrkred)  write(unitfile,format_smear) e(i), (sf_k(iz,0,i)), (sf_k(iz,1,i))
    enddo
 
   close(unit=unitfile)

  enddo 



!enddo



 deallocate(w_k)
 deallocate(sf_k)


end subroutine spectralfunction_r_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kvecs_hand (k_vecs2d, nk2d,nkx, nky)
real(kind=idp),allocatable, intent(out)	::k_vecs2d(:,:)
integer, intent(out)	::nk2d
integer, intent(in)	::nkx, nky
integer	::ik,ikn


if (nkx==2 .and. nky==2 .and. lkshifted_r)then 
nk2d=2

allocate(k_vecs2d(0:nk2d-1,4))

k_vecs2d=0
k_vecs2d(0,1)=0.25d0
k_vecs2d(0,2)=0.25d0
k_vecs2d(0,4)=2.0d0
k_vecs2d(1,1)=0.25d0
k_vecs2d(1,2)=-0.25d0
k_vecs2d(1,4)=2.0d0

else if (nkx==1 .and. nky==1 .and. .not.  lkshifted_r)then 
nk2d=1

allocate(k_vecs2d(0:nk2d-1,4))

k_vecs2d=0
k_vecs2d(0,1)=0.0d0
k_vecs2d(0,2)=0.0d0
k_vecs2d(0,4)=1.0d0



else 
print *, "set_kvecs_kr_hand not set"
stop
endif




!stop

end subroutine set_kvecs_hand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine part_ratio(e_r_tot,w_r_tot, lscfr)
real (kind=idp),intent(in)		:: e_r_tot(0:total_size-1,0:nk_r-1)
complex(KIND=idpc),intent(in)	 	:: w_r_tot(0:total_size-1,0:total_size-1,0:nk_r-1)
real (kind=idp), allocatable		:: pr(:,:)
real (kind=idp)				:: w2, w4, w2s
logical					::lscfr

allocate(pr(0:total_size-1,0:nk_r-1))

if (lscfr) OPEN(unit=16,file=trim(label)//'/e_r_scf',status='unknown')
if (.not. lscfr .and. lpathnscfr) OPEN(unit=16,file=trim(label)//'/e_r_path',status='unknown')
if (.not. lscfr .and. .not. lpathnscfr) OPEN(unit=16,file=trim(label)//'/e_r_nscf',status='unknown')


do ik=0,nk_r-1
 do i=0, total_size-1		!eigenvalue 
  w2=0
  w4=0
  do ix=0,nx-1	!x
   do iy=0,ny-1	!y
    do iz=0,nz-1	!z
     w2s=0
     do icf=0,1	!c/f
      do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     w2s=w2s+abs(w_r_tot(indp,i,ik))**2 !total on-site weight

      enddo
     enddo
     w2=w2+w2s
     w4=w4+w2s**2
    enddo
   enddo
  enddo
  if (abs(w2-1)>1e-5) print *, "part ratio w2 .ne. 1", i,ik,indp
pr(i,ik)=w2/w4
  write (16, "(i4,i7,2e15.6)") ik,i,e_r_tot(i,ik), pr(i,ik)
 enddo
enddo





deallocate(pr)
 close(unit=16)


end subroutine part_ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectralfunction_r_xyzk_zk(ik,lpathnscfr)
integer	,intent(in)	::ik
real(kind(ids)), allocatable	:: SF_k(:,:,:)
complex(kind=ids),allocatable:: rot_kx(:,:,:), wk(:,:)
 CHARACTER(2) :: ax, ay, az
 CHARACTER(3) :: ak, ak2
 CHARACTER(5) :: aka
 logical, intent(in)::lpathnscfr
real(kind(ids)), allocatable	:: SF_r(:,:,:,:,:)
real(kind(ids)), allocatable	:: SF_r2(:,:,:,:,:)
integer	::unitfile, ik2
real (kind=idp)	::weight_k
integer				::ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8
real(kind(ids)), allocatable	:: e_new(:), w_new(:,:), e_old(:), w_old(:,:)
integer	::n_states_new, n_states_old, nkx_l, nky_l
real(kind(ids)):: k_vec_arpes(2),k_vec_arpes_eq(2)
integer	::sx,sy,swap_xy, ysameasx

!if (ik==0) print *, n_k_total
!if (ik==0) print *, "sfk"


!!!!!!!
!PATH

if (lpathnscfr) then
 allocate(wk(0:4*nz_plot-1,0:total_size-1))
 allocate(sf_k(0:nz_plot-1,0:1,0:total_size-1))

 if (ik==0 .and. lfullbz)  write(*, "(3a5, 2a10, a5,4a10)") "iktot", "nK_x","nK_y","K_x","K_y", "ik", "k_x", "k_y", "k_tot x", "k_tot y"

 if (lfullbz) write(*, "(3i5, 2f10.5, i5, 4f10.5)")   ik, list_kk(ik,1),list_kk(ik,2),real(list_kk(ik,1))/nx,real(list_kk(ik,2))/ny,&
				list_kk(ik,3),k_vecs_r(list_kk(ik,3),1)/nx, k_vecs_r(list_kk(ik,3),2)/nx,&
				 k_tot(ik,1), k_tot(ik,2) 


 
!eigenvectors in kx,ky,z basis

wk=0

 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
   do i=0,total_size-1
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
!     print *, iz
     do icf=0,1	!c/f
      do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind2p=iz_ind + icf*nz_plot + is*2*nz_plot

if (.not. lwfctot)  wk(ind2p,i)= wk(ind2p,i)+w(indp,i)                    * exp(-2*(0,1)*pi*( k_tot(ik,1)*ix/real(1)+k_tot(ik,2)*iy/real(1) ) )/sqrt(real(nx*ny))
if (lwfctot)        wk(ind2p,i)= wk(ind2p,i)+w_r_tot(indp,i,list_kk(ik,3))* exp(-2*(0,1)*pi*( k_tot(ik,1)*ix/real(1)+k_tot(ik,2)*iy/real(1) ) )/sqrt(real(nx*ny))


       enddo 
      enddo
     enddo
    enddo
   enddo
  enddo


!spectral weights for k points

sf_k=0

   do i=0,total_size-1
    do icf=0,1	!c/f
     do is=0,1	!up/down 
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

!print *, indp, 4*nz_plot, i, total_size
      indp=iz_ind + icf*nz_plot + is*2*nz_plot
      sf_k(iz_ind,icf,i)=sf_k(iz_ind,icf,i)+real((conjg(wk(indp,i))*wk(indp,i)))		!sum over spin

      enddo
     enddo
    enddo
   enddo



!!!!!!
!write
!k resolved, average over z

    WRITE(UNIT=ak,FMT='(I3.3)') ik


    unitfile=100
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_k'//ak,status='unknown')

    write(unitfile,"(i10)") total_size

    do i=0, total_size-1		!eigenvalue 
if (lwfctot)          write(unitfile,format_smear) e_r_tot(i,list_kk(ik,3)), sum(sf_k(:,0,i))/real(nz), sum(sf_k(:,1,i))/real(nz)
if (.not. lwfctot)    write(unitfile,format_smear) e(i), sum(sf_k(:,0,i))/real(nz), sum(sf_k(:,1,i))/real(nz)
    enddo
    
   close(unit=unitfile)

!k-z resolved

     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)


!!!!!!!!!!!!!!
!COMPRESS

    n_states_old=total_size
    allocate(e_old(0:n_states_old-1))
    allocate(w_old(0:1,0:n_states_old-1))

    if (lwfctot) e_old(:)=e_r_tot(:,list_kk(ik,3))
    if (.not. lwfctot) e_old(:)=e(:)
    w_old(:,:)=sf_k(iz,:,:)

    call compress(e_old,w_old, n_states_old,2, e_new, w_new, n_states_new)
!!!!!!!!

    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    unitfile=100+iz
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nz'//az//'_k'//ak,status='unknown')

!not compressed
   ! write(unitfile,"(i10)") total_size
   !  do i=0, total_size-1		!eigenvalue 
   ! if (lwfctot)  write(unitfile,format_smear) e_r_tot(i,list_kk(ik,3)), (sf_k(iz,0,i)), (sf_k(iz,1,i))
   ! if (.not. lwfctot)  write(unitfile,format_smear) e(i), (sf_k(iz,0,i)), (sf_k(iz,1,i))
   ! enddo

!compressed
    write(unitfile,"(i10)") n_states_new
    do i=0, n_states_new-1		!eigenvalue 
     write(unitfile,format_smear) e_new(i), w_new(0,i), w_new(1,i)
    enddo


 deallocate(e_old)
 deallocate(w_old)
 deallocate(e_new)
 deallocate(w_new)


 
   close(unit=unitfile)

  enddo 


 deallocate(wk)
 deallocate(sf_k)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NSCF

else

allocate(SF_r(0:nx-1,0:ny-1,0:nz_plot-1, 0:2, 0:total_size-1))

!0:c
!1:f
!2:cf

SF_r=0
!weight_k=k_vecs_r(ik,4)/nk_r_tot
weight_k=1

 do i=0, total_size-1		!eigenvalue 
  do ix=0,nx-1	!x
   do iy=0,ny-1	!y
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
     do icf=0,1	!c/f
      do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind =ix + iy*nx + iz*nx*ny 

if (icf==0 .and. lwfctot)    sf_r(ix,iy,iz,icf,i)=sf_r(ix,iy,iz,icf,i)+abs(w_r_tot(indp,i,ik))**2*weight_k
if (icf==1 .and. lwfctot)    sf_r(ix,iy,iz,icf,i)=sf_r(ix,iy,iz,icf,i)+abs(w_r_tot(indp,i,ik))**2*weight_k
if (icf==0 .and. .not. lwfctot)    sf_r(ix,iy,iz,icf,i)=sf_r(ix,iy,iz,icf,i)+abs(w(indp,i))**2*weight_k
if (icf==1 .and. .not. lwfctot)    sf_r(ix,iy,iz,icf,i)=sf_r(ix,iy,iz,icf,i)+abs(w(indp,i))**2*weight_k

        icf2=1-icf
        do is2=0,1
         ind2p=ix + iy*nx + iz*nx*ny + icf2*nx*ny*nz + is2*2*nx*ny*nz
if (lwfctot)         sf_r(ix,iy,iz,2,i)=sf_r(ix,iy,iz,2,i)+real(conjg(w_r_tot(indp,i,ik))*(w_r_tot(ind2p,i,ik)))*weight_k
if (.not. lwfctot)   sf_r(ix,iy,iz,2,i)=sf_r(ix,iy,iz,2,i)+real(conjg(w(indp,i))*(w(ind2p,i)))*weight_k
        enddo
     
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
!!!!!!!!

!correction for reduced k set: no , to be done in postprocessing

!if (lkred_r .and. nx==ny) then

!allocate(SF_r2(0:nx-1,0:ny-1,0:nz-1, 0:1, 0:total_size-1))

! do ix=0,nx-1	!x
!  do iy=0,ny-1	!y

!sf_r2(ix,iy,:,:,:)=	(sf_r(ix,iy,:,:,:)+&
!			 sf_r(nx-1-ix,iy,:,:,:)+&
!			 sf_r(ix,ny-1-iy,:,:,:)+&
!			 sf_r(nx-1-ix,ny-1-iy,:,:,:)+&
!			 sf_r(iy,ix,:,:,:)+&
!			 sf_r(ny-1-iy,ix,:,:,:)+&
!			 sf_r(iy,nx-1-ix,:,:,:)+&
!			 sf_r(ny-1-iy,nx-1-ix,:,:,:))/8


!  enddo
! enddo

!sf_r=sf_r2

!deallocate(sf_r2)

!endif


!!!!!
!print

 do ix=0,nx-1	!x
  do iy=0,ny-1	!y
     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)
   
    n_states_old=total_size
    allocate(e_old(0:n_states_old-1))
    allocate(w_old(0:2,0:n_states_old-1))

    if (lwfctot) e_old(:)=e_r_tot(:,ik)
    if (.not. lwfctot) e_old(:)=e(:)
    w_old(:,:)=sf_r(ix,iy,iz,:,:)

!COMPRESS
    call compress(e_old,w_old, n_states_old,3, e_new, w_new, n_states_new)


    WRITE(UNIT=ax,FMT='(I2.2)') ix
    WRITE(UNIT=ay,FMT='(I2.2)') iy
    WRITE(UNIT=az,FMT='(I2.2)') iz
    WRITE(UNIT=ak,FMT='(I3.3)') ik
 
    ind=ix +iy*nx +iz*nx*ny
    unitfile=100+ind
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_nx'//ax//'_ny'//ay//'_nz'//az//'_k'//ak,status='unknown')

    write(unitfile,"(i10)") n_states_new
    do i=0, n_states_new-1		!eigenvalue 
     write(unitfile,format_smear) e_new(i), w_new(0,i), w_new(1,i), w_new(2,i)
    enddo
     write(unitfile,format_smear) b(ind)
     write(unitfile,format_smear) k_vecs_r(ik,1),k_vecs_r(ik,2),k_vecs_r(ik,4)/nk_r_tot

!not compressed
!    write(unitfile,"(i10)") total_size
!    do i=0, total_size-1		!eigenvalue 
!if (lwfctot)     write(unitfile,format_smear) e_r_tot(i,ik), sf_r(ix,iy,iz,0,i), sf_r(ix,iy,iz,1,i), sf_r(ix,iy,iz,2,i)
!if (.not. lwfctot)     write(unitfile,format_smear) e(i),    sf_r(ix,iy,iz,0,i), sf_r(ix,iy,iz,1,i), sf_r(ix,iy,iz,2,i)
!    enddo
!     write(unitfile,format_smear) b(ind)




 
 deallocate(e_old)
 deallocate(w_old)
 deallocate(e_new)
 deallocate(w_new)
 
    close(unit=unitfile)
 
   enddo
  enddo
 enddo


deallocate(sf_r)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NSCF, ARPES energy resolved

write (*,"(a5,3a3,2a5,6a10)" ), "ik", "sx","sy","sw", "ikx","iky","k_v_x","k_v_y","k_ar_x","k_ar_y","k_ar_x_eq","k_ar_y_eq"

if (mod(nkx_r_nscf,2)==1) then
do ik_arpes=0, nk_r_arpes-1
 if(list_k_vecs_r_arpes(ik_arpes,1)==ik) then
 
 ikx=list_k_vecs_r_arpes(ik_arpes,2)
 iky=list_k_vecs_r_arpes(ik_arpes,3)
 sx=list_k_vecs_r_arpes(ik_arpes,4)
 sy=list_k_vecs_r_arpes(ik_arpes,5)
 swap_xy=list_k_vecs_r_arpes(ik_arpes,6)
 k_vec_arpes(1)=k_vecs_r_arpes(ik_arpes,1)
 k_vec_arpes(2)=k_vecs_r_arpes(ik_arpes,2)
 k_vec_arpes_eq(1)=(k_vecs_r(ik,1)+ikx)/nx	!vectors to which calculation is applied
 k_vec_arpes_eq(2)=(k_vecs_r(ik,2)+iky)/nx



 allocate(wk(0:4*nz_plot-1,0:total_size-1))
 allocate(sf_k(0:nz_plot-1,0:1,0:total_size-1))


write (*,"(i5,3i3,2i5,6f10.5)" ), ik, sx,sy,swap_xy, ikx,iky,k_vecs_r(ik,1),k_vecs_r(ik,2),k_vec_arpes(1),k_vec_arpes(2),k_vec_arpes_eq(1),k_vec_arpes_eq(2)

      wk=0

     do ix=0,nx-1	!x
      do iy=0,ny-1	!y
       do i=0,total_size-1
        do iz_ind=0,nz_plot-1	!z
         iz=z_values_plot(iz_ind)
!     print *, iz
         do icf=0,1	!c/f
          do is=0,1	!up/down 

     indp=ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz
     ind2p=iz_ind + icf*nz_plot + is*2*nz_plot

  !    print *, "a",i,total_size,ind2p, 4*nz_plot


!if (.not. lwfctot)  wk(ind2p,i)= wk(ind2p,i)+w(indp,i)         * exp(-2*(0,1)*pi*( k_vec_arpes(1)*ix/real(1)+k_vec_arpes(2)*iy/real(1) ) )/sqrt(real(nx*ny))
!if (lwfctot)        wk(ind2p,i)= wk(ind2p,i)+w_r_tot(indp,i,ik)* exp(-2*(0,1)*pi*( k_vec_arpes(1)*ix/real(1)+k_vec_arpes(2)*iy/real(1) ) )/sqrt(real(nx*ny))
if (.not. lwfctot)  wk(ind2p,i)= wk(ind2p,i)+w(indp,i)         * exp(-2*(0,1)*pi*( k_vec_arpes_eq(1)*ix/real(1)+k_vec_arpes_eq(2)*iy/real(1) ) )/sqrt(real(nx*ny))
if (lwfctot)        wk(ind2p,i)= wk(ind2p,i)+w_r_tot(indp,i,ik)* exp(-2*(0,1)*pi*( k_vec_arpes_eq(1)*ix/real(1)+k_vec_arpes_eq(2)*iy/real(1) ) )/sqrt(real(nx*ny))


           enddo 
          enddo
         enddo
        enddo
       enddo
      enddo

sf_k=0


   do i=0,total_size-1
    do icf=0,1	!c/f
     do is=0,1	!up/down 
      do iz_ind=0,nz_plot-1	!z
      iz=z_values_plot(iz_ind)

      indp=iz_ind + icf*nz_plot + is*2*nz_plot
 !     print *, "b",i,total_size,indp, 4*nz_plot
      sf_k(iz_ind,icf,i)=sf_k(iz_ind,icf,i)+real((conjg(wk(indp,i))*wk(indp,i)))		!sum over spin

      enddo
     enddo
    enddo
   enddo
!stop
!!!!!!!!!!!!!!
!COMPRESS

     do iz_ind=0,nz_plot-1	!z
     iz=z_values_plot(iz_ind)

    n_states_old=total_size
    allocate(e_old(0:n_states_old-1))
    allocate(w_old(0:1,0:n_states_old-1))

    if (lwfctot) e_old(:)=e_r_tot(:,ik)
    if (.not. lwfctot) e_old(:)=e(:)
    w_old(:,:)=sf_k(iz,:,:)

    call compress(e_old,w_old, n_states_old,2, e_new, w_new, n_states_new)
!!!!!!!!
!writes

    WRITE(UNIT=az,FMT='(I2.2)') iz
    WRITE(UNIT=aka,FMT='(I5.5)') ik_arpes
 
    unitfile=100+iz
    OPEN(unit=unitfile,file=trim(label)//'/sf/sf_r_arpes_nz'//az//'_k'//aka,status='unknown')

!not compressed
   ! write(unitfile,"(i10)") total_size
   !  do i=0, total_size-1		!eigenvalue 
   ! if (lwfctot)  write(unitfile,format_smear) e_r_tot(i,list_kk(ik,3)), (sf_k(iz,0,i)), (sf_k(iz,1,i))
   ! if (.not. lwfctot)  write(unitfile,format_smear) e(i), (sf_k(iz,0,i)), (sf_k(iz,1,i))
   ! enddo

!compressed
    write(unitfile,"(i10)") n_states_new
    do i=0, n_states_new-1		!eigenvalue 
     write(unitfile,format_smear) e_new(i), w_new(0,i), w_new(1,i)
    enddo
     write(unitfile,format_smear) k_vec_arpes(1),k_vec_arpes(2), 1.0d0


 deallocate(e_old)
 deallocate(w_old)
 deallocate(e_new)
 deallocate(w_new)


 
   close(unit=unitfile)
enddo

 deallocate(wk)
 deallocate(sf_k)

endif

    enddo



endif

 




endif







end subroutine spectralfunction_r_xyzk_zk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compress(e, w, n_e,dim_w2, e_cpr, w_cpr, n_e_cpr)
real(kind(idp)),intent(in)	:: e(0:n_e-1)
integer, intent(in)	::n_e, dim_w2
real(kind(ids)),intent(in)	:: w(0:dim_w2-1,0:n_e-1)
real(kind(ids)),intent(out),allocatable	:: w_cpr(:,:)
real(kind(ids)),intent(out),allocatable	:: e_cpr(:)
integer, intent(out)	::n_e_cpr
integer	::c
real(kind(idp))	:: thr

thr=T/10

 c=0
i=0
do while (i<n_e-1) 
 j=0
 do while (i+j <n_e-1 .and. abs(e(i+j)-e(i))< thr) 
  j=j+1
 ! print *, i,j,c
 enddo
 i=i+j
 c=c+1
enddo

n_e_cpr=c


allocate(e_cpr(0:n_e_cpr-1))
allocate(w_cpr(0:dim_w2-1,0:n_e_cpr-1))

e_cpr=0
w_cpr=0

 c=0
i=0
do while (i<n_e-1) 
 j=0
 e_cpr(c)=e(i+j)!+thr/2
 do while (i+j <n_e-1 .and. abs(e(i+j)-e(i))< thr) 
  w_cpr(:,c)=w_cpr(:,c)+w(:,i+j)
  j=j+1
 enddo
 i=i+j
 c=c+1
enddo


!e_cpr=e
!w_cpr=W
!print *, "a"

end subroutine compress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_zvalues(nz_plot, z_values_plot)
integer, intent(inout)	::nz_plot
integer, allocatable, intent(inout)	:: z_values_plot(:)
integer	::iz

!if(nz_plot>5) nz_plot=nz
allocate (z_values_plot(0:nz_plot-1))

!if (nz_plot==nz) then
!do iz=0, nz-1
! z_values_plot(iz)=iz
!enddo
!else if (nz_plot==1) then
!z_values_plot(0)=0
!else if (nz_plot==2) then
!z_values_plot(0)=0
!z_values_plot(1)=1
!else if (nz_plot==3) then
!z_values_plot(0)=0
!z_values_plot(1)=nz-1
!z_values_plot(2)=nz/2
!else if (nz_plot==4) then
!z_values_plot(0)=0
!z_values_plot(1)=nz-1
!z_values_plot(2)=nz/2
!z_values_plot(3)=1
!else if (nz_plot==5) then
!z_values_plot(0)=0
!z_values_plot(1)=nz-1
!z_values_plot(2)=nz/2
!z_values_plot(3)=1
!z_values_plot(4)=nz-2
!else 
!do iz=0, nz-1
! z_values_plot(iz)=iz
!enddo
!endif
do iz=0, nz_plot-1
 z_values_plot(iz)=iz
enddo



print *, "z points for sf:", nz_plot
do iz=0, nz_plot-1
 print *,  z_values_plot(iz)
enddo


end subroutine set_zvalues


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine do_green
complex(kind=idpc),allocatable:: G0(:,:), G(:,:), V(:,:), Id(:,:),G0_1(:,:),G_1(:,:), G0V(:,:),  G0_x_tmp(:,:,:), G_born(:,:), V_kz(:,:), V_kknn(:,:), V_xz(:,:), V_knn(:,:,:)
complex(kind=idpc),allocatable:: G0_k(:,:,:), G0_x(:,:), dos_fourier(:,:,:),G0_x1(:,:),G0_x2(:,:), g_b(:,:), g0_b(:,:),G0_k_b(:,:,:),G0_k_d(:,:),Id_kr(:,:), Ham_kr(:,:), G0_1_b(:,:),&
G0_x_tmp_b(:,:,:),G0_x1_b(:,:),G0_x2_b(:,:), G0_x_b(:,:), Ham0_b(:,:), Ham_b(:,:), G_1_b(:,:), H_k_x(:,:), H_k_x_tmp(:,:,:),H_k_x_b(:,:), exp_factor(:,:),G0_k_e(:,:),G0_x_e(:,:),G_ew_k(:,:),&
exp_factor2(:,:), dos_r(:,:), dos_k(:,:),dos_ke(:,:),dos_ko(:,:),dos_re(:,:),dos_ro(:,:),exp_factore(:,:), exp_factoro(:,:),dos_fouriere(:,:,:),dos_fouriero(:,:,:), V_red(:,:), Id_red(:,:), dos_r_fft(:),&
g_fourier(:,:,:)
complex(kind=idpc),allocatable:: temp_1(:,:),temp_2(:,:),temp_3(:,:),temp_4(:,:),temp_5(:,:),temp_6(:,:), kn_xyzsa(:,:), kn_kzsa(:,:), kzsa_xyzsa(:,:), rhovrho(:,:),dos_born(:), dos_0(:),&
				 dos_g(:), dos_cf(:,:,:), dos_cf_born(:,:,:), dos_cf_0(:,:,:), dos_cf_r(:,:,:,:)
complex(kind=idpc),allocatable:: Ham(:,:), Ham0(:,:), T_m(:,:), temp(:,:), G0T(:,:), Ham0_gk(:,:),w(:,:),w0(:,:), G_wwa(:,:), G_wwb(:,:), G_ew(:,:), u(:,:), w_k_d(:,:), G_k_d(:,:), &
temp_k(:,:),G0_k_fft(:),G0_x_fft_tot(:,:,:), S_m(:,:), exphi(:),exphi_t(:),exphi_v(:), jdos(:),jdos2(:), &
jdos_klarge(:), ados_p(:),ados_p2(:), ados_m(:),simm_arpes(:,:,:,:,:),simm_arpes_klarge(:,:,:,:,:),&
g0_arpes(:,:,:,:)
real(kind=idp),allocatable:: b_kr0(:), lambda_kr0(:), dos_fin(:,:,:,:), b_total(:,:,:), lambda_total(:,:,:), e(:), e0(:),rho_kn(:,:,:), rhorho(:,:), dos_arpes(:,:,:,:), spin_arpes(:,:,:,:,:),&
				t_hop(:), g_fin(:,:,:,:), dos_arpes_klarge(:,:,:,:),spin_arpes_klarge(:,:,:,:,:),spin_arpes_rot(:,:,:,:,:)
real(kind=idp):: mu_kr0, v2, exp_corrx, exp_corry, xcenter, ycenter,exp_corrx_kr,exp_corry_kr, norm
integer	::indp1, indp2,indp3,indp4, size_large, ind3,ind1,ind4, unitfile, size_red, iz2_ind,  indptot,i2,iq,ixyz, nr_even, nr_odd, unit_write, ix_temp, &
indp5, indp6,indp7,indp8,nxy_red, size_v_red,indp9,indp10,ix3,iy3,ix4,iy4,ix5,iy5, iymax,ixmin1,ixmax1,iymin1,iymax1, ixcenter, iycenter, ixplotmin, ixplotmax, iyplotmin, iyplotmax,ixmax, index_kr,&
entry_dos, ind6p, ind7p
 character(7)	 ::ae2
 CHARACTER(2) :: ax, ay, az
 CHARACTER(4) :: ax4, ay4, az3
 logical:: ldiag_gg0, ltg0gg0, all_g, all_g_fourier
real(kind=idp):: time_0, time_1,time_gg0t, time_g0k, time_g0r_r, time_g0rr, time_g0tg0, time_dos, time_0_w,time_1_w, index_v(0:3,0:3,0:1), time_fft, time_fft_0, time_fft_1
complex(kind=idpc):: factor
integer :: status = 0, ignored_status
integer :: l_dfti(2),ind3b,ind3c,ind3d,ind3e,ind3f
!type(DFTI_DESCRIPTOR), POINTER :: hand


!hand => null()

!exp factors for corrections
if (mod(nkx_plot,2)==0) exp_corrx=-0.5d0
if (mod(nkx_plot,2)==1) exp_corrx=-0.5d0*dble(nkx_plot-1)/dble(nkx_plot)
if (mod(nky_plot,2)==0) exp_corry=-0.5d0
if (mod(nky_plot,2)==1) exp_corry=-0.5d0*dble(nky_plot-1)/dble(nky_plot)

if (mod(nkx_plot*nx_kr,2)==0) exp_corrx_kr=-0.5d0
if (mod(nkx_plot*nx_kr,2)==1) exp_corrx_kr=-0.5d0*dble(nkx_plot*nx_kr-1)/dble(nkx_plot*nx_kr)
if (mod(nky_plot*ny_kr,2)==0) exp_corry_kr=-0.5d0
if (mod(nky_plot*ny_kr,2)==1) exp_corry_kr=-0.5d0*dble(nky_plot*ny_kr-1)/dble(nky_plot*ny_kr)


!exp_corrx=0
!exp_corry=0
!exp_corrx_kr=0
!exp_corry_kr=0


lc4v=0
lc2v=0
lcs=0
lnos=1
!if(nx_kr==ny_kr) then
!lc4v=1
!lnos=0
!else
!lcs=1
!lnos=0
!endif

print *, "nx_kr=", nx_kr
print *, "ny_kr=", ny_kr
print *, "nx=", nx
print *, "ny=", ny
print *, "nky_plot=", nky_plot
print *, "nkx_plot=", nkx_plot
print *, "c_4v=", lc4v
print *, "c_2v=", lc2v
print *, "c_s=", lcs
print *, "no sym=", lnos


 call system('mkdir ' //trim(label)//"/g")

!ldiag_g0k=1	!diag rather than inversion for G0k large

time_gg0t=0
time_g0k=0
time_g0r_r=0
time_g0rr=0
time_g0tg0=0
time_dos=0
time_write=0
time_fourier=0
time_vkn=0
time_rhovrho=0
time_fft=0

if(lc4v .or. lc2v) then
if (mod(nkx_plot,2)==1) nx_red=nkx_plot/2+1
if (mod(nky_plot,2)==1) ny_red=nky_plot/2+1
if (mod(nkx_plot,2)==0) nx_red=nkx_plot/2
if (mod(nky_plot,2)==0) ny_red=nky_plot/2
elseif(lcs) then
if (mod(nkx_plot,2)==1) nx_red=nkx_plot
if (mod(nky_plot,2)==1) ny_red=nky_plot/2+1
if (mod(nkx_plot,2)==0) nx_red=nkx_plot
if (mod(nky_plot,2)==0) ny_red=nky_plot/2
elseif(lnos) then
if (mod(nkx_plot,2)==1) nx_red=nkx_plot
if (mod(nky_plot,2)==1) ny_red=nky_plot
if (mod(nkx_plot,2)==0) nx_red=nkx_plot
if (mod(nky_plot,2)==0) ny_red=nky_plot
endif

print *, "nx_red=", nx_red
print *, "ny_red=", ny_red

print *, "Doing Green"



allocate (t_hop(0:dim_hk-1))
if (l_orbital==1) then
t_hop(0)=t_d1
t_hop(1)=t_d1
t_hop(2)=-t_f1	!*b_k
t_hop(3)=t_f1	!*b_k
else if (l_orbital==2) then
t_hop(0)=t_d1
t_hop(1)=t_d1
t_hop(2)=t_d2
t_hop(3)=t_d2
t_hop(4)=-t_f1	!*b_k
t_hop(5)=t_f1	!*b_k
t_hop(6)=-t_f2	!*b_k
t_hop(7)=t_f2	!*b_k
endif
if (l_orbital_f==3) then
t_hop(8)=-t_f3	!*b_k
t_hop(9)=t_f3	!*b_k
endif

write (*,*) "STM t parameters:"
do i=0, dim_hk-1
write (*, "(i5,10f10.5)"), i, t_hop(i)
enddo


print *, "H,H0,G,G_0,V"

!calculation only at Gamma
deallocate(k_vecs_r)
nk_r=1
allocate(k_vecs_r(0:nk_r-1, 4))
k_vecs_r(0,1)=0
k_vecs_r(0,2)=0
k_vecs_r(0,3)=0
k_vecs_r(0,4)=1

!H,H0,V

allocate(b_kr0(0:n_sites-1))
allocate(lambda_kr0(0:n_sites-1))
!allocate(G0(0:total_size-1,0:total_size-1))
allocate(Id  (0:total_size-1,0:total_size-1))
allocate(V   (0:total_size-1,0:total_size-1))
allocate(Ham (0:total_size-1,0:total_size-1))
allocate(Ham0(0:total_size-1,0:total_size-1))
allocate(G_ew(0:total_size-1,0:total_size-1))

size_v_red=2*nx*ny*nz_scf*(l_orbital+l_orbital_f)	!only the first nz_scf layers
allocate(V_red   (0:size_v_red-1,0:size_v_red-1))
allocate(Id_red  (0:size_v_red-1,0:size_v_red-1))


!identity
Id=0
do i=0, total_size-1
 Id(i,i)=1.0d0
enddo

Id_red=0
do i=0, size_v_red-1
 Id_red(i,i)=1.0d0
enddo

!H
do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
   ind=index_xyz(ix,iy,iz,nx,ny)
   index_kr=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr)
    b_kr0(ind)=b_kr(index_kr)
    lambda_kr0(ind)=lambda_kr(index_kr)
   enddo
  enddo
 enddo
mu_kr0=mu_kr



Ham=0
if (.not. lv_nscf .and. (lgreen .or. lborn)) call build_ham_r(0,Ham,b,lambda,mu,ec_site,ef_site)
if (lv_nscf  .and. (lgreen .or. lborn))      call build_ham_r(0,Ham,b_kr0,lambda_kr0,mu,ec_site,ef_site)
!H0
!here I remove the on-site energies different from the std value
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1

ind=index_xyz(ix,iy,iz,nx,ny)
index_kr=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr)
ef_site(ind,:)=ef_site_kr(index_kr,:)
ec_site(ind,:)=ec_site_kr(index_kr,:)

enddo
enddo
enddo

!print *, ef_site
!print *, ec_site
!print * , d_ec, d_ef, d_ef7
!print *, ef_site
!print *, ec_site

 call build_ham_r(0,Ham0,b_kr0,lambda_kr0,mu_kr0,ec_site,ef_site)


if (abs(mu-mu_kr0)>1e-3) then
print *, "different mu", mu, mu_kr0
endif


!V=H-H0
V=Ham-Ham0

!V_red( only in the first nz_scf layers)
V_red=0
 do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz_scf-1
    do icf=0,1
     do is=0,1
      do o1=0, l_orb(icf)-1
      ind= index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz)
      indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz_scf)
 do ix2=0,nx-1
  do iy2=0,ny-1
   do iz2=0,nz_scf-1
    do icf2=0,1
     do is2=0,1
      do o2=0, l_orb(icf2)-1
      ind2= index_xyzcso(ix2,iy2,iz2,icf2,is2,o2,nx, ny, nz)
      ind2p=index_xyzcso(ix2,iy2,iz2,icf2,is2,o2,nx, ny, nz_scf)
      
      V_red(indp, ind2p)=V(ind,ind2)

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


!for the 2d eff model: impurity prop to identity in spin space
if (l2dmodel) then
V_red=0
V_red(0,0)=d_ef
V_red(1,1)=V_red(0,0)
endif

!print *, d_ef
!print *, V_red

!writes V to be plotted (only diagonal part)
OPEN(unit=100,file=trim(label)//'/v_g',status='unknown')

do iz=0,nz-1	!z
 do ix=0,nx-1	!x
  do iy=0,ny-1	!y

      indp1=index_xyzcso(ix,iy,iz,0,0,0,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz	!c1 down
      indp2=index_xyzcso(ix,iy,iz,0,1,0,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz	!c1 up
      indp3=index_xyzcso(ix,iy,iz,0,0,1,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz	!c2 down
      indp4=index_xyzcso(ix,iy,iz,0,1,1,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz	!c2 up
      indp5=index_xyzcso(ix,iy,iz,1,0,0,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz + 4*nx*ny*nz	!f1 down
      indp6=index_xyzcso(ix,iy,iz,1,1,0,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz + 4*nx*ny*nz	!f1 up
      indp7=index_xyzcso(ix,iy,iz,1,0,1,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz + 4*nx*ny*nz	!f2 down
      indp8=index_xyzcso(ix,iy,iz,1,1,1,nx, ny, nz)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz + 4*nx*ny*nz	!f2 up
      indp9=index_xyzcso(ix,iy,iz,1,0,2,nx, ny, nz)	!f3 down
     indp10=index_xyzcso(ix,iy,iz,1,1,2,nx, ny, nz)	!f3 up

if (l_orbital==2 .and. l_orbital_f==2) write(100,"(3i5,16e10.3)"), ix,iy,iz, real(V(indp1,indp1)),real(V(indp2,indp2)),real(V(indp3,indp3)),real(V(indp4,indp4)),&
				     real(V(indp5,indp5)),real(V(indp6,indp6)),real(V(indp7,indp7)),real(V(indp8,indp8))!,&
!				     aimag(V(indp1,indp1)),aimag(V(indp2,indp2)),aimag(V(indp3,indp3)),aimag(V(indp4,indp4)),&
!				     aimag(V(indp5,indp5)),aimag(V(indp6,indp6)),aimag(V(indp7,indp7)),aimag(V(indp8,indp8))

if (l_orbital==2 .and. l_orbital_f==3) write(100,"(3i5,20e10.3)"), ix,iy,iz, real(V(indp1,indp1)),real(V(indp2,indp2)),real(V(indp3,indp3)),real(V(indp4,indp4)),&
				     real(V(indp5,indp5)),real(V(indp6,indp6)),real(V(indp7,indp7)),real(V(indp8,indp8)),real(V(indp9,indp9)),real(V(indp10,indp10))!,&
!				     aimag(V(indp1,indp1)),aimag(V(indp2,indp2)),aimag(V(indp3,indp3)),aimag(V(indp4,indp4)),&
!				     aimag(V(indp5,indp5)),aimag(V(indp6,indp6)),aimag(V(indp7,indp7)),aimag(V(indp8,indp8)),aimag(V(indp9,indp9)),aimag(V(indp10,indp10))

if (l_orbital==1 .and. l_orbital_f==1) write(100,"(3i5,20e10.3)"), ix,iy,iz, real(V(indp1,indp1)),real(V(indp2,indp2)),real(V(indp5,indp5)),real(V(indp6,indp6))!,&
!				     real(V(indp5,indp5)),real(V(indp6,indp6)),real(V(indp7,indp7)),real(V(indp8,indp8)),real(V(indp9,indp9)),real(V(indp10,indp10))!,&

!if (l_orbital==1) write(100,"(3i5,16e15.6)"), ix,iy,iz, real(V(indp1,indp1)),real(V(indp2,indp2)),real(V(indp3,indp3)),real(V(indp4,indp4)),&
!				     aimag(V(indp1,indp1)),aimag(V(indp2,indp2)),aimag(V(indp5,indp5)),aimag(V(indp6,indp6))
		

  enddo
 enddo
write(100,*), ""
enddo

!reduced
 close(unit=100)
OPEN(unit=100,file=trim(label)//'/v_g_red',status='unknown')
print *, "nz_scf", nz_scf
do iz=0,nz_scf-1	!z
 do ix=0,nx-1	!x
  do iy=0,ny-1	!y

      indp1=index_xyzcso(ix,iy,iz,0,0,0,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz	!c1 down
      indp2=index_xyzcso(ix,iy,iz,0,1,0,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz	!c1 up
      indp3=index_xyzcso(ix,iy,iz,0,0,1,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz	!c2 down
      indp4=index_xyzcso(ix,iy,iz,0,1,1,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz	!c2 up
      indp5=index_xyzcso(ix,iy,iz,1,0,0,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 0*2*nx*ny*nz + 4*nx*ny*nz	!f1 down
      indp6=index_xyzcso(ix,iy,iz,1,1,0,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 0*nx*ny*nz + 1*2*nx*ny*nz + 4*nx*ny*nz	!f1 up
      indp7=index_xyzcso(ix,iy,iz,1,0,1,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz + 4*nx*ny*nz	!f2 down
      indp8=index_xyzcso(ix,iy,iz,1,1,1,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz + 4*nx*ny*nz	!f2 up
      indp9=index_xyzcso(ix,iy,iz,1,0,2,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 0*2*nx*ny*nz + 4*nx*ny*nz	!f3 down
      indp10=index_xyzcso(ix,iy,iz,1,1,2,nx, ny, nz_scf)	!ix + iy*nx + iz*nx*ny + 1*nx*ny*nz + 1*2*nx*ny*nz + 4*nx*ny*nz	!f3 up

if (l_orbital==2  .and. l_orbital_f==3) write(100,"(3i5,20e15.6)"), ix,iy,iz, real(V_red(indp1,indp1)),real(V_red(indp2,indp2)),real(V_red(indp3,indp3)),real(V_red(indp4,indp4)),&
				     real(V_red(indp5,indp5)),real(V_red(indp6,indp6)),real(V_red(indp7,indp7)),real(V_red(indp8,indp8)),real(V_red(indp9,indp9)),real(V_red(indp10,indp10))!,&
!				     aimag(V_red(indp5,indp5)),aimag(V_red(indp6,indp6)),aimag(V_red(indp7,indp7)),aimag(V_red(indp8,indp8))
!				     aimag(V_red(indp1,indp1)),aimag(V_red(indp2,indp2)),aimag(V_red(indp3,indp3)),aimag(V_red(indp4,indp4)),&
if (l_orbital==2  .and. l_orbital_f==2) write(100,"(3i5,16e15.6)"), ix,iy,iz, real(V_red(indp1,indp1)),real(V_red(indp2,indp2)),real(V_red(indp3,indp3)),real(V_red(indp4,indp4)),&
				     real(V_red(indp5,indp5)),real(V_red(indp6,indp6)),real(V_red(indp7,indp7)),real(V_red(indp8,indp8))!,&
!				     aimag(V_red(indp5,indp5)),aimag(V_red(indp6,indp6)),aimag(V_red(indp7,indp7)),aimag(V_red(indp8,indp8))
!				     aimag(V_red(indp1,indp1)),aimag(V_red(indp2,indp2)),aimag(V_red(indp3,indp3)),aimag(V_red(indp4,indp4)),&
if (l_orbital==1) write(100,"(3i5,16e15.6)"), ix,iy,iz, real(V_red(indp1,indp1)),real(V_red(indp2,indp2)),real(V_red(indp5,indp5)),real(V_red(indp6,indp6)),&
				     aimag(V_red(indp1,indp1)),aimag(V_red(indp2,indp2)),aimag(V_red(indp5,indp5)),aimag(V_red(indp6,indp6))
if (l_orbital==0) write(100,"(3i5,16e15.6)"), ix,iy,iz, real(V_red(indp5,indp5)),real(V_red(indp6,indp6)),&
						aimag(V_red(indp5,indp5)),aimag(V_red(indp6,indp6))


  enddo
 enddo
write(100,*), ""
enddo

 close(unit=100)

!stop

!deallocate(Id)
!deallocate(Ham)
!deallocate(Ham0)
deallocate(b_kr0)
deallocate(lambda_kr0)



!!!!!!!!loop over omega

size_large=nkx_plot*nky_plot*nz*2*(l_orbital+l_orbital_f)*nx_kr*ny_kr
if (lc4v) nxy_red=(nx_red*(nx_red+1))/2	!C_4v: just in one octant (y equiv to x)
if (lc2v.or.lcs.or.lnos) nxy_red=nx_red*ny_red	!C_2v: just in one quadrant !C_s: just in the lower half; nx_red=nkx_plot!no sym: everywhere
!if (lcs)  nxy_red=nx_red*ny_red	
!if (lnos)  nxy_red=nx_red*ny_red
size_red=nxy_red*nz_plot*2*(l_orbital+l_orbital_f)*nx_kr*ny_kr

print *, "nxy_red=",nxy_red

print *, "nz_scf:", nz_scf
print *, "nz_plot:", nz_plot
print *, "nz_g0:", nz_g0
print *, "Total_size(small):", total_size
print *, "Size_large:", size_large
print *, "Nx_red:", nx_red
print *, "Nx_red*(Nx_red+1)/2:", (nx_red*(nx_red+1))/2
print *, "Size_red:", size_red		!by C4v symmetry, and nz_plot=1 (where I finally compute G)
print *, "Size_v_red:", size_v_red	!by nz_scf


do omega=omega_min,omega_max,delta_omega

all_g=1
all_g_fourier=0
if (l2dmodel) all_g=0
!if(abs(omega)<1e-5) all_g=1
!if (lqpi) all_g=1

WRITE(UNIT=ae2,FMT='(SPf7.4)') omega

lqpi=1
if (omega < en_min_qpi-0.0001 .or. omega > en_max_qpi+0.0001) lqpi=0

write (*, "(a10,f10.5)") "omega=", omega
print *, "qpi=", lqpi
if (lv_nscf) print *, "WARNING:NSCF"



!!!!!!!!!!!!!!!!!!!!!!!!G0_k(saz,s'a'z',k)
print *, "G0k"

time_0=secnds(0.0)

allocate(G0_k(0:size_h_kr_red-1,0:size_h_kr_red-1,0:nkr_plot-1))
G0_k=0

!if(ldiag_g0k) then
!with eigenvectors

allocate(G_ew_k(0:size_h_kr_red_e-1,0:size_h_kr_red_e-1))
allocate(temp_k(0:size_h_kr_red-1,0:size_h_kr_red_e-1))
allocate(w_k_d(0:size_h_kr_red-1,0:size_h_kr_red_e-1))	!eigevectors:indp (zsa), i (eigenvalue)
allocate(G_k_d(0:size_h_kr_red-1,0:size_h_kr_red-1))

G0_k=0
 do ik=0, nkr_plot-1
 
  G_ew_k=0
  do i=0, size_h_kr_red_e-1
   G_ew_k(i,i)=1/(omega-e_kr_tot_nscf(i,ik)+(0,1)*T_g)
  enddo
  
 temp_k=0
 G_k_d=0
 w_k_d(:,:)=w_kr_tot_nscf(:,:,ik)

! call matmul1(size_h_kr_red, w_k_d, G_ew_k, temp_k)
! call matmul_hc(size_h_kr_red, temp_k, w_k_d, G_k_d)
 
 !w_k_d(size_h_kr_red,size_h_kr) * G_ew_k(size_h_kr,size_h_kr) = temp_k (size_h_kr_red, size_h_kr)
 call matmul2(size_h_kr_red,size_h_kr_red_e,size_h_kr_red_e, w_k_d, G_ew_k, temp_k)

 !temp_k(size_h_kr_red,size_h_kr) * w_k_d hc(size_h_kr,size_h_kr_red) = G_k_d (size_h_kr_red, size_h_kr_red)
 call matmul2hc(size_h_kr_red,size_h_kr_red,size_h_kr_red_e, temp_k, w_k_d, G_k_d)
 G0_k(:,:,ik)=G_k_d(:,:)
 
enddo

deallocate(G_ew_k)
deallocate(temp_k)
deallocate(w_k_d)
deallocate(G_k_d)


!!!! writes G0_k on file -->arpes in kx-ky plane (z=0)
allocate(dos_arpes(0:nkr_plot-1,0:1, 0:nz_plot-1, 0:l_orbital_f-1)) !real
allocate(g0_arpes(0:nkr_plot-1,0:1, 0:nz_plot-1, 0:l_orbital_f-1))	!complex

dos_arpes=0
g0_arpes=0

do ik=0, nkr_plot-1
do ix=0,nx_kr-1
do iy=0,ny_kr-1
  do iz=0,nz_plot-1			!layer
   do icf=0,1			!c/f
    do is=0,1			!up/down sum) 
     do o1=0,l_orb(icf)-1
     
!   indp=index_zcso(iz,icf,is,o1,nz_scf)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
   indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_plot)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
   indp2=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf

dos_arpes(ik,icf,iz,o1)=dos_arpes(ik,icf,iz,o1)-1/pi*aimag(G0_k(indp2,indp2,ik))
g0_arpes (ik,icf,iz,o1)=g0_arpes (ik,icf,iz,o1)+(G0_k(indp2,indp2,ik))

enddo
enddo
enddo
enddo
enddo
enddo
enddo


if (lqpi .and. nx_kr==1 .and. ny_kr==1) then
do iz=0,min(3,nz_plot-1)	!z
 WRITE(UNIT=az,FMT='(I2.2)') iz
 OPEN(unit=100,file=trim(label)//'/g/a_kr_arpesxy_nz'//az//'_e'//ae2,status='unknown')

do ik=0, nkr_plot-1
!c1, c2, f1, f2
if (l_orbital==2 .and. l_orbital_f==2) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), dos_arpes(ik,0,iz,0), dos_arpes(ik,0,iz,1),dos_arpes(ik,1,iz,0),dos_arpes(ik,1,iz,1)
!c1, c2, f1, f2, f3
if (l_orbital==2 .and. l_orbital_f==3) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), &
					dos_arpes(ik,0,iz,0), dos_arpes(ik,0,iz,1),dos_arpes(ik,1,iz,0),dos_arpes(ik,1,iz,1),dos_arpes(ik,1,iz,2)
!c,f
if (l_orbital==1) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), dos_arpes(ik,0,iz,0) ,dos_arpes(ik,1,iz,0)

!only f (2dmodel)
if (l_orbital==0) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), dos_arpes(ik,1,iz,0)

enddo !ik

 close(unit=100)
enddo	!iz

!if (nz_scf>1) then
 OPEN(unit=100,file=trim(label)//'/g/a_kr_arpesxy_nztot_e'//ae2,status='unknown')

do ik=0, nkr_plot-1
!c1, c2, f1, f2
if (l_orbital==2 .and. l_orbital_f==2) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), sum(dos_arpes(ik,0,:,0)), sum(dos_arpes(ik,0,:,1)),&
				     									     sum(dos_arpes(ik,1,:,0)), sum(dos_arpes(ik,1,:,1))
!c1, c2, f1, f2, f3
if (l_orbital==2 .and. l_orbital_f==3) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), &
					sum(dos_arpes(ik,0,:,0)), sum(dos_arpes(ik,0,:,1)),sum(dos_arpes(ik,1,:,0)),sum(dos_arpes(ik,1,:,1)),sum(dos_arpes(ik,1,:,2))
!c,f
if (l_orbital==1) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), sum(dos_arpes(ik,0,:,0)) ,sum(dos_arpes(ik,1,:,0))

!only f (2dmodel)
if (l_orbital==0) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), sum(dos_arpes(ik,1,:,0))

enddo	!ik
!endif	!nz_scf

 close(unit=100)


endif


!JDOS
if (ljdos) then
print *, "JDOS"
allocate(ados_p(0:nkx_plot*nky_plot-1))
allocate(ados_p2(0:nkx_plot*nky_plot-1))
!allocate(ados_m(0:nkx_plot*nky_plot-1))
allocate(jdos(0:nkx_plot*nky_plot-1))
allocate(jdos2(0:nkx_plot*nky_plot-1))

ados_p=0
do ik=0, nkr_plot-1
  do iz=0,nz_plot-1			!layer
   do icf=0,1			!c/f
     do o1=0,l_orb(icf)-1
ados_p(ik)=ados_p(ik)+dos_arpes(ik,icf,iz,o1)		!for simplicity I sum over orbitals and layers
ados_p2(ik)=ados_p2(ik)+g0_arpes(ik,icf,iz,o1)		!for simplicity I sum over orbitals and layers
!ados_p2(ik)=aimag(ados_p2(ik))
!ados_m(ik)=ados_m(ik)+dos_arpes(ik,icf,iz,o1)
enddo
enddo
enddo
enddo

norm=0
norm2=0
do ik=0, nkr_plot-1
norm=norm+abs(ados_p(ik))**2
norm2=norm2+abs(ados_p2(ik))**2
enddo

!do ikx=0,nkx_plot-1	
!do iky=0,nky_plot-1	
!iq=ikx+iky*nkx_plot
!ados_p(iq)=ados_p(iq)/nky_plot/nkx_plot!*exp(-2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!ados_m(iq)=ados_m(iq)/nky_plot/nkx_plot*exp(+2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!enddo
!enddo


l_dfti(1)=nkx_plot
l_dfti(2)=nky_plot

 ! status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

 ! status = DftiCommitDescriptor(hand)
 !if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, ados_p)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"
!  status = DftiComputeForward(hand, ados_p2)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"


jdos=0
jdos2=0
do ikx=0,nkx_plot-1	
do iky=0,nky_plot-1	
iq=ikx+iky*nkx_plot
!ados_p(iq)=ados_p(iq)/nky_plot/nkx_plot*exp(-2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!ados_m(iq)=ados_m(iq)/nky_plot/nkx_plot*exp(+2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
jdos(iq)=abs(ados_p(iq))**2*exp(-2*pi*(0,1)*(ikx*exp_corrx+iky*exp_corry))/nkx_plot/nky_plot!corr fct to center around (n-1)/2n in x space !*ados_m(ix)		
jdos2(iq)=abs(ados_p2(iq))**2*exp(-2*pi*(0,1)*(ikx*exp_corrx+iky*exp_corry))/nkx_plot/nky_plot!corr fct to center around (n-1)/2n in x space !*ados_m(ix)		
enddo
enddo



l_dfti(1)=nkx_plot
l_dfti(2)=nky_plot

 ! status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
 !if (0 /= status) print *, "error"

 ! status = DftiCommitDescriptor(hand)
 !if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, jdos)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"
!  status = DftiComputeForward(hand, jdos2)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"

 ! status = DftiFreeDescriptor(hand)
 !if (0 /= status) print *, "error"


!do ikx=0,nkx_plot-1	
!do iky=0,nky_plot-1	
!iq=ikx+iky*nkx_plot
!jdos(iq)=jdos(iq)*exp(-2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!enddo
!enddo


OPEN(unit=100,file=trim(label)//'/g/jdos_e'//ae2,status='unknown')
do iq=0, nkx_plot*nky_plot-1	!now in real space
write(100,"(2f10.5,16e16.8)"), kr_vecs_nscf(iq,5),kr_vecs_nscf(iq,6), abs(jdos(iq))/norm,real(jdos2(iq))/norm2,aimag(jdos2(iq))/norm2!, real(jdos(iq)),aimag(jdos(iq))
enddo
write(100,"(16e16.8)"), norm,norm2
 close(unit=100)


deallocate(jdos)
deallocate(ados_p)
deallocate(jdos2)
deallocate(ados_p2)
!deallocate(ados_m)
endif !ljdos


!ARPES in the extended BZ !only works for 2x1 reconstruction
if (nx_kr>1 .or. ny_kr>1) then
allocate(dos_arpes_klarge(0:nkr_plot_klarge-1,0:1, 0:nz_plot-1, 0:l_orbital_f-1))	
dos_arpes_klarge=0

do ik=0, nkr_plot_klarge-1
do ix=0,nx_kr-1
do iy=0,ny_kr-1
do ix2=0,nx_kr-1
do iy2=0,ny_kr-1
  do iz=0,nz_plot-1			!layer
   do icf=0,1			!c/f
    do is=0,1			!up/down sum) 
     do o1=0,l_orb(icf)-1
     
   indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
   indp2=index_xyzcso(ix2,iy2,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf

dos_arpes_klarge(ik,icf,iz,o1)=dos_arpes_klarge(ik,icf,iz,o1)-1/pi*aimag(G0_k(indp,indp2,list_kk(ik,3))*&
			exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo

if (lqpi) then
do iz=0,min(3,nz_plot-1)	!z
 WRITE(UNIT=az,FMT='(I2.2)') iz
 OPEN(unit=100,file=trim(label)//'/g/a_kr_arpesek_nz'//az//'_e'//ae2,status='unknown')

do ik=0, nkr_plot_klarge-1
!c1, c2, f1, f2
if (l_orbital==2 .and. l_orbital_f==2) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), dos_arpes_klarge(ik,0,iz,0),&
								dos_arpes_klarge(ik,0,iz,1),dos_arpes_klarge(ik,1,iz,0),dos_arpes_klarge(ik,1,iz,1)
!c1, c2, f1, f2, f3
if (l_orbital==2 .and. l_orbital_f==3) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), &
					dos_arpes_klarge(ik,0,iz,0), dos_arpes_klarge(ik,0,iz,1),dos_arpes_klarge(ik,1,iz,0),dos_arpes_klarge(ik,1,iz,1),dos_arpes_klarge(ik,1,iz,2)
!c,f
if (l_orbital==1) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), dos_arpes_klarge(ik,0,iz,0) ,dos_arpes_klarge(ik,1,iz,0)

!only f (2dmodel)
if (l_orbital==0) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), dos_arpes_klarge(ik,1,iz,0)

enddo

 close(unit=100)
enddo


!if(nz_scf>1) then
 OPEN(unit=100,file=trim(label)//'/g/a_kr_arpesek_nztot_e'//ae2,status='unknown')

do ik=0, nkr_plot_klarge-1
!c1, c2, f1, f2
if (l_orbital==2 .and. l_orbital_f==2) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), sum(dos_arpes_klarge(ik,0,:,0)),&
								sum(dos_arpes_klarge(ik,0,:,1)),sum(dos_arpes_klarge(ik,1,:,0)),sum(dos_arpes_klarge(ik,1,:,1))
!c1, c2, f1, f2, f3
if (l_orbital==2 .and. l_orbital_f==3) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), &
					sum(dos_arpes_klarge(ik,0,:,0)), sum(dos_arpes_klarge(ik,0,:,1)),sum(dos_arpes_klarge(ik,1,:,0)),sum(dos_arpes_klarge(ik,1,:,1)),sum(dos_arpes_klarge(ik,1,:,2))
!c,f
if (l_orbital==1) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), sum(dos_arpes_klarge(ik,0,:,0)) ,sum(dos_arpes_klarge(ik,1,:,0))

!only f (2dmodel)
if (l_orbital==0) write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), sum(dos_arpes_klarge(ik,1,:,0))

enddo	!k
!endif	!nz_scf

 close(unit=100)

endif



!JDOS
if (ljdos) then
print *, "JDOS_rec"
allocate(ados_p(0:nkr_plot_klarge-1))
!allocate(ados_m(0:nkx_plot*nky_plot-1))
allocate(jdos(0:nkr_plot_klarge-1))

ados_p=0
do ik=0, nkr_plot_klarge-1
  do iz=0,nz_plot-1			!layer
   do icf=0,1			!c/f
     do o1=0,l_orb(icf)-1
ados_p(ik)=ados_p(ik)+dos_arpes_klarge(ik,icf,iz,o1)		!for simplicity I sum over orbitals and layers
enddo
enddo
enddo
enddo

 OPEN(unit=100,file=trim(label)//'/g/a_kr_arpesek_b_nztot_e'//ae2,status='unknown')
do ik=0, nkr_plot_klarge-1
if (l_orbital==0) write(100,"(2f10.5,16e17.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), abs(ados_p(ik))
enddo
 close(unit=100)
 
norm=0
do ik=0,nkr_plot_klarge-1
norm=norm+abs(ados_p(ik))**2
enddo

l_dfti(1)=nkx_plot*nx_kr
l_dfti(2)=nky_plot*ny_kr

!  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

!  status = DftiCommitDescriptor(hand)
! if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, ados_p)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"

jdos=0
do ikx=0,nkx_plot*nx_kr-1	
do iky=0,nky_plot*ny_kr-1	
iq=ikx+iky*nkx_plot*nx_kr
!ados_p(iq)=ados_p(iq)/nky_plot/nkx_plot*exp(-2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!ados_m(iq)=ados_m(iq)/nky_plot/nkx_plot*exp(+2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
jdos(iq)=abs(ados_p(iq))**2/nkr_plot_klarge*exp(-2*pi*(0,1)*(ikx*exp_corrx_kr+iky*exp_corry_kr))!corr fct to center around (n-1)/2n in x space !*ados_m(ix)		
enddo
enddo



!  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

!  status = DftiCommitDescriptor(hand)
! if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, jdos)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"


!do ikx=0,nkx_plot-1	
!do iky=0,nky_plot-1	
!iq=ikx+iky*nkx_plot
!jdos(iq)=jdos(iq)*exp(-2*pi*(0,1)*(exp_corrx*ikx+exp_corry*iky))	!additional phase factor from having k points from -0.5 instead from 0
!enddo
!enddo


OPEN(unit=100,file=trim(label)//'/g/jdos_ek_e'//ae2,status='unknown')
do iq=0, nkr_plot_klarge-1	
write(100,"(2f10.5,16e16.8)"), kr_vecs_nscf_klarge(iq,1),kr_vecs_nscf_klarge(iq,2), abs(jdos(iq))/norm!, real(jdos(iq)),aimag(jdos(iq))
enddo
 close(unit=100)

deallocate(jdos)
deallocate(ados_p)
endif !jdos


deallocate(dos_arpes_klarge)


endif 	!extended BZ


!!!!!!!!!!!!!spin structure
if (lqpi .and. nkx_plot*nx_kr <100) then
!sum over orbital degeneracy, c/f resolved
print *, "spin structure"

allocate(spin_arpes_rot(0:nkr_plot-1,0:1,1:3, 0:nz_plot-1,0:3))
allocate(spin_arpes(0:nkr_plot-1,0:1,1:3, 0:nz_plot-1,0:3))
allocate(simm_arpes(0:nkr_plot-1,0:0,0:n_sym_arpes,0:0,0:0))
!			k        c/f xyz     z layer  components
spin_arpes=0
simm_arpes=0
spin_arpes_rot=0

do ik=0, nkr_plot-1
 do iz=0,nz_plot-1			!layer
do ix=0, nx_kr-1
do iy=0, ny_kr-1
  do icf=0,1			!c/f
   do is=0,1			!up/down (sum) 
    do is2=0,1
      do o1=0, l_orb(icf)-1	!?	!OK, for d if 2 it gives 0
       do o2=0, l_orb(icf)-1
     
   indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
   indp2=index_xyzcso(ix,iy,iz,icf,is2,o2,nx_kr,ny_kr,nz_g0)   !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
     
!   indp=index_zcso(iz,icf,is,o1,nz_scf)! iz + icf*nz_scf + is*2*nz_scf  + o1*4*nz_scf
!   indp2=index_zcso(iz,icf,is2,o2,nz_scf)!iz + icf*nz_scf + is2*2*nz_scf + o2*4*nz_scf

     do ixyz=1,3
spin_arpes(ik,icf,ixyz,iz,0)=spin_arpes(ik,icf,ixyz,iz,0)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) 
if (icf==1 .and. o1==2 .and. o2==2) spin_arpes(ik,icf,ixyz,iz,1)=spin_arpes(ik,icf,ixyz,iz,1)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) !contribution from Gamma7 only
if (icf==1 .and. o1<2 .and. o2<2) spin_arpes(ik,icf,ixyz,iz,2)=spin_arpes(ik,icf,ixyz,iz,2)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) !contribution from Gamma8 only
if (icf==1 .and. ((o1<2 .and. o2==2) .or. (o1==2 .and. o2<2))) spin_arpes(ik,icf,ixyz,iz,3)=spin_arpes(ik,icf,ixyz,iz,3)-1/pi*aimag((G0_k(indp,indp2,ik))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) !contribution from Gamma7-Gamma8 mixed
enddo


do ixyz=0,n_sym_arpes
simm_arpes(ik,0,ixyz,0,0)=simm_arpes(ik,0,ixyz,0,0)-1/pi*(G0_k(indp,indp2,ik)*simm_refl(is2,is,o2,o1,icf,ixyz)) !no aimag
enddo

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


do ik=0, nkr_plot-1
 do iz=0,nz_plot-1			!layer
  do icf=0,1			!c/f
do ixyz=1,3
do ixyz2=1,3
spin_arpes_rot(ik,icf,ixyz,iz,0)=spin_arpes_rot(ik,icf,ixyz,iz,0)+rot_mat_surf(ixyz,ixyz2)*spin_arpes(ik,icf,ixyz2,iz,0)
enddo
enddo
enddo
enddo
enddo



!WRITE(UNIT=ae2,FMT='(SPf7.4)') omega

!do iz=0,min(1,nz_scf-1)
do iz=0,nz_plot-1
 WRITE(UNIT=az,FMT='(I2.2)') iz
 OPEN(unit=100,file=trim(label)//'/g/a_kr_spins_nz'//az//'_e'//ae2,status='unknown')
! OPEN(unit=101,file=trim(label)//'/g/a_kr_spins_rot_nz'//az//'_e'//ae2,status='unknown')

do ik=0, nkr_plot-1

if(l2dmodel) then
write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), 	spin_arpes(ik,1,1,iz,0), spin_arpes(ik,1,2,iz,0), spin_arpes(ik,1,3,iz,0),&	!fx,fy,fz, 345
									sum(dos_arpes(ik,1,iz,:))							!dos f
else

!write(100,"(2f10.5,30f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), 	&
!write(100,"(2f10.5,30f15.8)"), k_surf(1,1)*kr_vecs_nscf(ik,1)+k_surf(2,1)*kr_vecs_nscf(ik,2),k_surf(1,2)*kr_vecs_nscf(ik,1)+k_surf(2,2)*kr_vecs_nscf(ik,2),&
!									spin_arpes(ik,0,1,iz,0), spin_arpes(ik,0,2,iz,0), spin_arpes(ik,0,3,iz,0),&	!cx,cy,cz, 345
!									spin_arpes(ik,1,1,iz,0), spin_arpes(ik,1,2,iz,0), spin_arpes(ik,1,3,iz,0),&	!fx,fy,fz, 678
!									spin_arpes(ik,1,1,iz,1), spin_arpes(ik,1,2,iz,1), spin_arpes(ik,1,3,iz,1),&	!f7x,f7y,f7z
!									spin_arpes(ik,1,1,iz,2), spin_arpes(ik,1,2,iz,2), spin_arpes(ik,1,3,iz,2),&	!f8x,f8y,f8z
!									spin_arpes(ik,1,1,iz,3), spin_arpes(ik,1,2,iz,3), spin_arpes(ik,1,3,iz,3),&	!f78x,f78y,f78z
!									sum(dos_arpes(ik,0,iz,:)),sum(dos_arpes(ik,1,iz,:))				!dos c, dos f

write(100,"(2f10.5,30f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), &!k_surf(1,1)*kr_vecs_nscf(ik,1)+k_surf(2,1)*kr_vecs_nscf(ik,2),k_surf(1,2)*kr_vecs_nscf(ik,1)+k_surf(2,2)*kr_vecs_nscf(ik,2),&
!write(101,"(2f10.5,30f15.8)"), k_surf(1,1)*kr_vecs_nscf(ik,1)+k_surf(2,1)*kr_vecs_nscf(ik,2),k_surf(1,2)*kr_vecs_nscf(ik,1)+k_surf(2,2)*kr_vecs_nscf(ik,2),&
									spin_arpes_rot(ik,0,1,iz,0), spin_arpes_rot(ik,0,2,iz,0), spin_arpes_rot(ik,0,3,iz,0),&	!cx,cy,cz, 345
									spin_arpes_rot(ik,1,1,iz,0), spin_arpes_rot(ik,1,2,iz,0), spin_arpes_rot(ik,1,3,iz,0),&	!fx,fy,fz, 678
									spin_arpes_rot(ik,1,1,iz,1), spin_arpes_rot(ik,1,2,iz,1), spin_arpes_rot(ik,1,3,iz,1),&	!f7x,f7y,f7z
									spin_arpes_rot(ik,1,1,iz,2), spin_arpes_rot(ik,1,2,iz,2), spin_arpes_rot(ik,1,3,iz,2),&	!f8x,f8y,f8z
									spin_arpes_rot(ik,1,1,iz,3), spin_arpes_rot(ik,1,2,iz,3), spin_arpes_rot(ik,1,3,iz,3),&	!f78x,f78y,f78z
									sum(dos_arpes(ik,0,iz,:)),sum(dos_arpes(ik,1,iz,:))				!dos c, dos f


endif

enddo	!k

 close(unit=100)
! close(unit=101)
enddo	!z


!summed over z
 OPEN(unit=100,file=trim(label)//'/g/a_kr_spins_nztot_e'//ae2,status='unknown')

do ik=0, nkr_plot-1

write(100,"(2f10.5,30f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), 	sum(spin_arpes_rot(ik,0,1,:,0)), sum(spin_arpes_rot(ik,0,2,:,0)), sum(spin_arpes_rot(ik,0,3,:,0)),&	!cx,cy,cz, 345
									sum(spin_arpes_rot(ik,1,1,:,0)), sum(spin_arpes_rot(ik,1,2,:,0)), sum(spin_arpes_rot(ik,1,3,:,0)),&	!fx,fy,fz, 678
									sum(spin_arpes_rot(ik,1,1,:,1)), sum(spin_arpes_rot(ik,1,2,:,1)), sum(spin_arpes_rot(ik,1,3,:,1)),&	!f7x,f7y,f7z
									sum(spin_arpes_rot(ik,1,1,:,2)), sum(spin_arpes_rot(ik,1,2,:,2)), sum(spin_arpes_rot(ik,1,3,:,2)),&	!f8x,f8y,f8z
									sum(spin_arpes_rot(ik,1,1,:,3)), sum(spin_arpes_rot(ik,1,2,:,3)), sum(spin_arpes_rot(ik,1,3,:,3)),&	!f8x,f8y,f8z
									sum(dos_arpes(ik,0,:,:)),sum(dos_arpes(ik,1,:,:))				!dos c, dos f
enddo

 close(unit=100)


!!simm
 OPEN(unit=100,file=trim(label)//'/g/a_kr_simms_e'//ae2,status='unknown')

do ik=0, nkr_plot-1

!if(l2dmodel) then
!write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), 	spin_arpes(ik,1,1,iz,0), spin_arpes(ik,1,2,iz,0), spin_arpes(ik,1,3,iz,0),&	!fx,fy,fz, 345
!									sum(dos_arpes(ik,1,iz,:))							!dos f
!else


!eigenvalues of M are proportional to i, so i take the real part!-->I remove I from the definition: imaginary part
write(100,"(2f10.5,60f15.8)"), kr_vecs_nscf(ik,5),kr_vecs_nscf(ik,6), (aimag(simm_arpes(ik,0,j,0,0))/aimag(simm_arpes(ik,0,0,0,0)), j=1,n_sym_arpes),aimag(simm_arpes(ik,0,0,0,0))	

enddo	!k

 close(unit=100)






deallocate(spin_arpes)
deallocate(spin_arpes_rot)
deallocate(simm_arpes)

endif


!spin in the extended BZ
if (lqpi .and. nkx_plot <90 .and. (nx_kr>1 .or. ny_kr>1)) then
!sum over orbital degeneracy, c/f resolved

allocate(spin_arpes_klarge(0:nkr_plot_klarge-1,0:1,1:3, 0:nz_plot-1,0:3))
allocate(simm_arpes_klarge(0:nkr_plot_klarge-1,0:0,0:n_sym_arpes,0:0,0:0))
!			k        c/f xyz     z layer  components
spin_arpes_klarge=0

do ik=0, nkr_plot_klarge-1
 do iz=0,nz_plot-1			!layer
do ix=0, nx_kr-1
do iy=0, ny_kr-1
do ix2=0, nx_kr-1
do iy2=0, ny_kr-1
  do icf=0,1			!c/f
   do is=0,1			!up/down (sum) 
    do is2=0,1
      do o1=0, l_orb(icf)-1	!OK, for d if 2 it gives 0
       do o2=0, l_orb(icf)-1
     
   indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)   
   indp2=index_xyzcso(ix2,iy2,iz,icf,is2,o2,nx_kr,ny_kr,nz_g0)   
     
!   indp=index_zcso(iz,icf,is,o1,nz_scf)! iz + icf*nz_scf + is*2*nz_scf  + o1*4*nz_scf
!   indp2=index_zcso(iz,icf,is2,o2,nz_scf)!iz + icf*nz_scf + is2*2*nz_scf + o2*4*nz_scf

     do ixyz=1,3

spin_arpes_klarge(ik,icf,ixyz,iz,0)=spin_arpes_klarge(ik,icf,ixyz,iz,0)-1/pi*aimag(((G0_k(indp,indp2,list_kk(ik,3)))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz)) *&
							exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))

if (icf==1 .and. o1==2 .and. o2==2) spin_arpes_klarge(ik,icf,ixyz,iz,1)=spin_arpes_klarge(ik,icf,ixyz,iz,1)-1/pi*aimag(((G0_k(indp,indp2,list_kk(ik,3)))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz))*&
			exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))
 !contribution from Gamma7 only
if (icf==1 .and. o1<2 .and. o2<2) spin_arpes_klarge(ik,icf,ixyz,iz,2)=spin_arpes_klarge(ik,icf,ixyz,iz,2)-1/pi*aimag(((G0_k(indp,indp2,list_kk(ik,3)))*sigma_xyz_cf(is2,is,o2,o1,icf, ixyz))*&
			exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))
 !contribution from Gamma8 only
if (icf==1 .and. ((o1<2 .and. o2==2) .or. (o1==2 .and. o2<2))) spin_arpes_klarge(ik,icf,ixyz,iz,3)=spin_arpes_klarge(ik,icf,ixyz,iz,3)-1/pi*aimag(((G0_k(indp,indp2,list_kk(ik,3)))*&
sigma_xyz_cf(is2,is,o2,o1,icf, ixyz))*exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))
 !contribution from Gamma7-Gamma8 mixed

enddo

do ixyz=0,n_sym_arpes
simm_arpes_klarge(ik,0,ixyz,0,0)=simm_arpes_klarge(ik,0,ixyz,0,0)-1/pi*(G0_k(indp,indp2,list_kk(ik,3))*simm_refl(is2,is,o2,o1,icf,ixyz)*&
							exp(2*pi*(0,1)*(kr_vecs_nscf_klarge(ik,1)*(ix-ix2)/nx_kr+kr_vecs_nscf_klarge(ik,2)*(iy-iy2)/ny_kr)))
enddo


enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo

!WRITE(UNIT=ae2,FMT='(SPf7.4)') omega

!do iz=0,min(1,nz_scf-1)
do iz=0,nz_plot-1
 WRITE(UNIT=az,FMT='(I2.2)') iz
 OPEN(unit=100,file=trim(label)//'/g/a_kr_spins_ek_nz'//az//'_e'//ae2,status='unknown')

do ik=0, nkr_plot_klarge-1

if(l2dmodel) then
write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), spin_arpes_klarge(ik,1,1,iz,0), spin_arpes_klarge(ik,1,2,iz,0), spin_arpes_klarge(ik,1,3,iz,0)
else

write(100,"(2f10.5,30f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), spin_arpes_klarge(ik,0,1,iz,0), spin_arpes_klarge(ik,0,2,iz,0), spin_arpes_klarge(ik,0,3,iz,0),&	!cx,cy,cz, 345
									spin_arpes_klarge(ik,1,1,iz,0), spin_arpes_klarge(ik,1,2,iz,0), spin_arpes_klarge(ik,1,3,iz,0),&	!fx,fy,fz, 678
									spin_arpes_klarge(ik,1,1,iz,1), spin_arpes_klarge(ik,1,2,iz,1), spin_arpes_klarge(ik,1,3,iz,1),&	!f7x,f7y,f7z
									spin_arpes_klarge(ik,1,1,iz,2), spin_arpes_klarge(ik,1,2,iz,2), spin_arpes_klarge(ik,1,3,iz,2),&	!f8x,f8y,f8z
									spin_arpes_klarge(ik,1,1,iz,3), spin_arpes_klarge(ik,1,2,iz,3), spin_arpes_klarge(ik,1,3,iz,3)	
endif

enddo	!ik

 close(unit=100)
enddo	!iz

 OPEN(unit=100,file=trim(label)//'/g/a_kr_spins_ek_nztot_e'//ae2,status='unknown')

do ik=0, nkr_plot_klarge-1

write(100,"(2f10.5,30f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), 	&
									sum(spin_arpes_klarge(ik,0,1,:,0)), sum(spin_arpes_klarge(ik,0,2,:,0)), sum(spin_arpes_klarge(ik,0,3,:,0)),&	!cx,cy,cz, 345
									sum(spin_arpes_klarge(ik,1,1,:,0)), sum(spin_arpes_klarge(ik,1,2,:,0)), sum(spin_arpes_klarge(ik,1,3,:,0)),&	!fx,fy,fz, 678
									sum(spin_arpes_klarge(ik,1,1,:,1)), sum(spin_arpes_klarge(ik,1,2,:,1)), sum(spin_arpes_klarge(ik,1,3,:,1)),&	!f7x,f7y,f7z
									sum(spin_arpes_klarge(ik,1,1,:,2)), sum(spin_arpes_klarge(ik,1,2,:,2)), sum(spin_arpes_klarge(ik,1,3,:,2)),&	!f8x,f8y,f8z
									sum(spin_arpes_klarge(ik,1,1,:,3)), sum(spin_arpes_klarge(ik,1,2,:,3)), sum(spin_arpes_klarge(ik,1,3,:,3))	!f8x,f8y,f8z
enddo

 close(unit=100)


!!simm
 OPEN(unit=100,file=trim(label)//'/g/a_kr_simms_ek_e'//ae2,status='unknown')

do ik=0, nkr_plot_klarge-1

!if(l2dmodel) then
!write(100,"(2f10.5,16f15.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), 	spin_arpes(ik,1,1,iz,0), spin_arpes(ik,1,2,iz,0), spin_arpes(ik,1,3,iz,0),&	!fx,fy,fz, 345
!									sum(dos_arpes(ik,1,iz,:))							!dos f
!else


!eigenvalues of M are proportional to i, so i take the real part!-->I remove I from the definition: imaginary part
write(100,"(2f10.5,60f15.8)"), kr_vecs_nscf_klarge(ik,1),kr_vecs_nscf_klarge(ik,2), (aimag(simm_arpes_klarge(ik,0,j,0,0))/aimag(simm_arpes_klarge(ik,0,0,0,0)), j=1,n_sym_arpes),&
				aimag(simm_arpes_klarge(ik,0,0,0,0))	

enddo	!k

 close(unit=100)


deallocate(spin_arpes_klarge)
deallocate(simm_arpes_klarge)

endif








deallocate(dos_arpes)
deallocate(g0_arpes)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


time_1=secnds(0.0)
time_g0k=time_g0k+time_1-time_0


!G0(r-r',z,z') FFT of G0k
print *, "G0(r-r',z,z')"

time_0=secnds(0.0)

!allocate(G0_k(0:size_h_kr-1,0:size_h_kr-1,0:nkr_plot-1))
allocate(G0_x_tmp(0:size_h_kr_red-1,0:size_h_kr_red-1,0:nkr_plot-1)) !G_0(r-r',z,z')
G0_x_tmp=0



time_fft_0=secnds(0.0)

allocate(G0_k_fft(0:nkx_plot*nky_plot-1))
G0_k_fft=0


do iz=0, nz_g0-1
do ix=0, nx_kr-1
do iy=0, ny_kr-1
 do icf=0,1
  do is=0,1
   do o1=0,l_orb(icf)-1
   indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx_kr,ny_kr,nz_g0)  
   
ind1=indp
do iz2=0, nz_g0-1
do ix2=0, nx_kr-1
do iy2=0, ny_kr-1
 do icf2=0,1
  do is2=0,1
   do o2=0,l_orb(icf2)-1
   indp2=index_xyzcso(ix2,iy2,iz2,icf2,is2,o2,nx_kr,ny_kr,nz_g0)   
!   indp2= index_zcso(iz2,icf2,is2,o2,nz_scf)!iz2 + icf2*nz_scf + is2*2*nz_scf + o2*4*nz_scf
ind2=indp2
!write(*, "(6i4)") iz,icf,is,iz2,icf2,is2

!do ikx=0,nkx_plot-1
! do iky=0,nky_plot-1

!   	ik=ikx+(nkx_plot)*iky
	
!if (mod(nkx_plot,2)==0)   kr_vecs_nscf (ik,1)=-0.5d0+real(ikx)/real(nkx_plot)
!if (mod(nkx_plot,2)==0)   kr_vecs_nscf (ik,2)=-0.5d0+real(iky)/real(nkx_plot)

!if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,1)=-0.5d0+real(ikx+0.5d0)/real(nkx_plot)
!if (mod(nkx_plot,2)==1)   kr_vecs_nscf (ik,2)=-0.5d0+real(iky+0.5d0)/real(nkx_plot)

G0_k_fft=0
G0_k_fft(:)=G0_k(ind1,ind2,:)		!G_0(k) for fixed xyzsax'y'z's'a'

!   enddo
!  enddo



l_dfti(1)=nkx_plot
l_dfti(2)=nky_plot

!  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

!  status = DftiCommitDescriptor(hand)
! if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, G0_k_fft)	!-		r->k eikr 	k->r e-ikr
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"

!G0_x_fft_tot(:,indp,indp2)=G0_k_fft(:)	!G_0(r,zsa,z's'a')

!if (mod(nkx_plot,2)==0) exp_corrx=-0.5d0
!if (mod(nkx_plot,2)==1) exp_corrx=-0.5d0*dble(nkx_plot-1)/dble(nkx_plot)
!if (mod(nky_plot,2)==0) exp_corry=-0.5d0
!if (mod(nky_plot,2)==1) exp_corry=-0.5d0*dble(nky_plot-1)/dble(nky_plot)
!exp_corrx=0.5d0
!exp_corry=0


do ix3=0, nkx_plot-1
 do iy3=0, nky_plot-1
!ind=ix3+iy3*nkx_plot
ind=index_g0xtmp(ix3,iy3,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)

!G_0(r,zsa,z's'a')
G0_x_tmp(indp,indp2,ind)=G0_k_fft(ind)/nky_plot/nkx_plot*exp(-2*pi*(0,1)*(exp_corrx*ix3+exp_corry*iy3))	!additional phase factor from having k points from -0.5 instead from 0

enddo
enddo

!999 print '("  Error, status = ",I0)', status



enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo

deallocate(G0_k_fft)
!deallocate(G0_x_fft_tot)

time_g0r_r=time_g0r_r+time_1-time_0



deallocate (G0_k)


time_fft_1=secnds(0.0)
time_fft=time_fft+time_fft_1-time_fft_0


!G0(r,r',z,z')

print *, "G0(r,r',z,z')"

time_0=secnds(0.0)

allocate(dos_0 (0:size_red-1))			!dos_0(z,r) surface
allocate(dos_cf_0(0:nxy_red*nz_plot*nx_kr*ny_kr-1,0:dim_hk-1,0:dim_hk-1))
allocate(dos_cf_r(0:nxy_red*nz_plot*nx_kr*ny_kr-1,0:dim_hk-1,0:dim_hk-1,5))
allocate(G0_x1(0:size_red-1,0:size_v_red-1))	!G_0(z,z',r,r') surface - scatt
allocate(G0_x2(0:size_v_red-1,0:size_red-1))	!G_0(z,z',r,r') scatt - surface
allocate(G0(0:size_v_red-1,0:size_v_red-1))	!G_0(z,z',r,r') scatt - scatt
!size_red=nxy_red*nz_plot*2*(l_orbital+l_orbital_f)*nx_kr*ny_kr
!size_v_red=2*nx*ny*nz_scf*(l_orbital+l_orbital_f)	!only the first nz_scf layers



dos_0=0
G0_x1=0
G0_x2=0
G0=0
dos_cf_0=0
dos_cf_r=0

 call set_ixmax(ixmax)	!nx_red-1
 do ix=0,ixmax	!index of kr cell nx_red=nkx_plot(/2)(+1) nkx_plot=nkxy_plot/nx_kr
  call set_iymax(iymax,ix)
  do iy=0,iymax
   do ix2=0,nx_kr-1	!index of site inside kr cell
    do iy2=0,ny_kr-1
     do iz_ind=0,nz_plot-1	!z
      iz=z_values_plot(iz_ind)
      do icf=0,1
       do is=0,1
        do o1=0,l_orb(icf)-1
      ind=index_xyzcso_red(ix,iy,ix2,iy2,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
!        write(*, "(20i5)"), ix,iy,ix2,iy2,iz_ind,icf,is,o1,ind
!        write(*, "(20i5)") ind
!      indp= index_zcso(iz,icf,is,o1, nz_scf) !iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf 
      indp= index_xyzcso(ix2,iy2,iz,icf,is,o1,nx_kr,ny_kr,nz_g0) !index of G0_x_tmp


ind3=index_g0xtmp(0,0,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)
ind3b=index_g0xtmp(1,0,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)
ind3c=index_g0xtmp(0,1,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)
ind3d=index_g0xtmp(1,1,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)
ind3e=index_g0xtmp(2,0,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)
ind3f=index_g0xtmp(0,2,nkx_plot,nky_plot)		!r-r'=0; index of x_plot, y_plot (cell index)

	
dos_0(ind)=G0_x_tmp(indp,indp,ind3)		!surface - surface

     do is2=0,1
       do icf2=0,1
      do o2=0,l_orb(icf2)-1
      ind2p= index_xyzcso(ix2,iy2,iz,icf2,is2,o2,nx_kr,ny_kr,nz_g0) !index of G0
 
      ind5p=index_xyz_red(ix,iy,ix2,iy2,iz_ind,nz_plot, nx_red, nxy_red)	!position of atom
      ind3p=index_cso(icf,is,o1) !icf+2*is+4*o1
      ind4p=index_cso(icf2,is2,o2) !icf2+2*is2+4*o2

dos_cf_0(ind5p,ind3p,ind4p)=dos_cf_0(ind5p,ind3p,ind4p)+ G0_x_tmp(indp,ind2p,ind3)	!same position, different orbitals
dos_cf_r(ind5p,ind3p,ind4p,1)=dos_cf_r(ind5p,ind3p,ind4p,1)+ G0_x_tmp(indp,ind2p,ind3b)	!different position, different orbitals
dos_cf_r(ind5p,ind3p,ind4p,2)=dos_cf_r(ind5p,ind3p,ind4p,2)+ G0_x_tmp(indp,ind2p,ind3c)	!different position, different orbitals
dos_cf_r(ind5p,ind3p,ind4p,3)=dos_cf_r(ind5p,ind3p,ind4p,3)+ G0_x_tmp(indp,ind2p,ind3d)	!different position, different orbitals
dos_cf_r(ind5p,ind3p,ind4p,4)=dos_cf_r(ind5p,ind3p,ind4p,4)+ G0_x_tmp(indp,ind2p,ind3e)	!different position, different orbitals
dos_cf_r(ind5p,ind3p,ind4p,5)=dos_cf_r(ind5p,ind3p,ind4p,5)+ G0_x_tmp(indp,ind2p,ind3f)	!different position, different orbitals

     enddo
    enddo
   enddo

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
!enddo
!enddo
!enddo
!enddo
!enddo

print *, "dos_0"
do iz_ind=0,nz_plot-1	!z
 iz=z_values_plot(iz_ind)
   do ix2=0,nx_kr-1	!but dependent from atom in the cell
    do iy2=0,ny_kr-1

WRITE(UNIT=az,FMT='(I2.2)') iz
WRITE(UNIT=ax,FMT='(I2.2)') ix2
WRITE(UNIT=ay,FMT='(I2.2)') iy2

    
OPEN(unit=100,file=trim(label)//'/g/g0_xkr'//ax//'_ykr'//ay//'_z'//az//'_e'//ae2,status='unknown')
!OPEN(unit=301,file=trim(label)//'/g/g0_x01y00_z'//az,status='unknown')
!OPEN(unit=302,file=trim(label)//'/g/g0_x00y01_z'//az,status='unknown')
!OPEN(unit=303,file=trim(label)//'/g/g0_x01y01_z'//az,status='unknown')
!OPEN(unit=304,file=trim(label)//'/g/g0_x02y00_z'//az,status='unknown')
!OPEN(unit=305,file=trim(label)//'/g/g0_x00y02_z'//az,status='unknown')
!OPEN(unit=401,file=trim(label)//'/g0_x01y00_z'//az,status='unknown')
!OPEN(unit=402,file=trim(label)//'/g0_x00y01_z'//az,status='unknown')
!OPEN(unit=403,file=trim(label)//'/g0_x01y01_z'//az,status='unknown')
!OPEN(unit=404,file=trim(label)//'/g0_x02y00_z'//az,status='unknown')
!OPEN(unit=405,file=trim(label)//'/g0_x00y02_z'//az,status='unknown')


!if (l_orbital==2 .and. l_orbital_f==3) then
!write (*, "(10(1x,a16))") "c1-", "c1+", "c2-", "c2+", "f1-", "f1+", "f2-", "f2+", "f3-", "f3+"
!elseif (l_orbital==1 .and. l_orbital_f==1) then
!write (*, "(10(1x,a16))") "c-", "c+", "f-", "f+"
!endif

!write (*, "(8(1x,a16))") "c1-", "c1+", "f1-", "f1+"
ix=ix2	! dos_0 is independent from position
iy=iy2 
 ind5p=index_xyz_red(ix,iy,ix2,iy2,iz_ind,nz_plot, nx_red, nxy_red)	!position of atom
      !do icf=0,1
      ! do is=0,1
      !  do o1=0,l_orb(icf)-1
      !do icf2=0,1
      ! do is2=0,1
      !  do o2=0,l_orb(icf)-1
      !ind=index_cso (icf,is,o1)
      !ind2=index_cso(icf2,is2,o2)
do ind=0,dim_hk-1
!write (*  , "(10(1x,f9.4,f9.4))", advance='no') (dos_cf_0(ind5p, ind,ind2))
write (100, "(10(1x,f12.4,f12.4))") (dos_cf_0(ind5p, ind,ind2),ind2=0,dim_hk-1)
!write (101, "(10(1x,f12.4,f12.4))") (dos_cf_r(ind5p, ind,ind2,1),ind2=0,dim_hk-1)
!write (102, "(10(1x,f12.4,f12.4))") (dos_cf_r(ind5p, ind,ind2,2),ind2=0,dim_hk-1)
!write (103, "(10(1x,f12.4,f12.4))") (dos_cf_r(ind5p, ind,ind2,3),ind2=0,dim_hk-1)
!write (104, "(10(1x,f12.4,f12.4))") (dos_cf_r(ind5p, ind,ind2,4),ind2=0,dim_hk-1)
!write (105, "(10(1x,f12.4,f12.4))") (dos_cf_r(ind5p, ind,ind2,5),ind2=0,dim_hk-1)
!enddo
!enddo
enddo
!write (100, "(a1)") ""
!enddo
!enddo
!enddo
!write (*  , "(a1)") ""
!print *, ""

OPEN(unit=200,file=trim(label)//'/g/g0_xkr'//ax//'_ykr'//ay//'_z'//az,status='unknown',position='append')
OPEN(unit=201,file=trim(label)//'/g0_xkr'//ax//'_ykr'//ay//'_z'//az,status='unknown',position='append')
OPEN(unit=301,file=trim(label)//'/g/g0_x01y00_z'//az,status='unknown',position='append')
  OPEN(unit=401,file=trim(label)//'/g0_x01y00_z'//az,status='unknown',position='append')
OPEN(unit=302,file=trim(label)//'/g/g0_x00y01_z'//az,status='unknown',position='append')
  OPEN(unit=402,file=trim(label)//'/g0_x00y01_z'//az,status='unknown',position='append')
OPEN(unit=303,file=trim(label)//'/g/g0_x01y01_z'//az,status='unknown',position='append')
  OPEN(unit=403,file=trim(label)//'/g0_x01y01_z'//az,status='unknown',position='append')
OPEN(unit=304,file=trim(label)//'/g/g0_x02y00_z'//az,status='unknown',position='append')
  OPEN(unit=404,file=trim(label)//'/g0_x02y00_z'//az,status='unknown',position='append')
OPEN(unit=305,file=trim(label)//'/g/g0_x00y02_z'//az,status='unknown',position='append')
  OPEN(unit=405,file=trim(label)//'/g0_x00y02_z'//az,status='unknown',position='append')



write (200, "(f15.5,300f12.4)") omega, (dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk), i=0, dim_hk**2-1)
write (201, "(f15.5,300f12.4)") omega, (dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk), i=0, dim_hk**2-1)
write (301, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,1), i=0, dim_hk**2-1)
write (401, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,1), i=0, dim_hk**2-1)
write (302, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,2), i=0, dim_hk**2-1)
write (402, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,2), i=0, dim_hk**2-1)
write (303, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,3), i=0, dim_hk**2-1)
write (403, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,3), i=0, dim_hk**2-1)
write (304, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,4), i=0, dim_hk**2-1)
write (404, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,4), i=0, dim_hk**2-1)
write (305, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,5), i=0, dim_hk**2-1)
write (405, "(f15.5,300f12.4)") omega, (dos_cf_r(ind5p,mod(i,dim_hk),i/dim_hk,5), i=0, dim_hk**2-1)




!   do ix2=0,nx_kr-1	!index of site inside kr cell
!    do iy2=0,ny_kr-1
!      do icf=0,1
!       do is=0,1
!        do o1=0,l_orb(icf)-1
!      ind=index_xyzcso_red(ix,iy,ix2,iy2,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
!write (*, "(10(1x,f8.4,f8.4))", advance='no') (dos_0(ind))
!enddo
!enddo
!enddo
!enddo
!enddo
!print *, ""


!do i=0,size_h_kr_red-1
!write (*, "(100(1x,f8.4,f8.4))") (g0_x_tmp(i,j,0), j=0,size_h_kr_red-1)
!enddo



!write (*, "(10(1x,f8.4,f8.4))") (dos_0(j), j=0,dim_hk-1)


!write (*, "(10(1x,f8.4,f8.4))") (dos_cf_0(0, ind,ind2), j=0,dim_hk-1)
!write (*, "(10(1x,f8.4,f8.4))") (dos_0(j), j=0,dim_hk-1)


!OPEN(unit=100,file=trim(label)//'/g/g0_e'//ae2,status='unknown')
!do i=0,dim_hk-1
!write (100, "(20f15.8)") (dos_cf_0(0, i,j), j=0,dim_hk-1)
!enddo
 close(unit=100)
 close(unit=200)
 close(unit=201)
 close(unit=301)
 close(unit=302)
 close(unit=303)
 close(unit=304)
 close(unit=305)
 close(unit=401)
 close(unit=402)
 close(unit=403)
 close(unit=404)
 close(unit=405)

enddo
enddo
enddo


!G0 with different indices for matrix-matrix multiplication

if (lgreen .or. lborn) then

!ixmin1=(nkx_plot-nx/nx_kr)/2			!(5-1)/2=2	index of cell
!ixmax1=(nkx_plot-nx/nx_kr)/2+nx/nx_kr-1		!2
!iymin1=(nky_plot-ny/ny_kr)/2			!(9-1)/2=4
!iymax1=(nky_plot-ny/ny_kr)/2+ny/ny_kr-1		!4
ixmin1=(nkx_plot*nx_kr-nx)/2			!20 4		index of site (total)
ixmax1=(nkx_plot*nx_kr-nx)/2+nx-1		!21 5
iymin1=(nky_plot*ny_kr-ny)/2			!20 (9-1)/2=4
iymax1=(nky_plot*ny_kr-ny)/2+ny-1		!20 4

!print *, "G0(r,r',z,z') left"

!G0_x left
!all x and y (one 1/8 of the surface); z only surface
 call set_ixmax(ixmax)	!ixmax=nx_red-1
 do ix=0,ixmax		!index of kr cell nx_red=nkx_plot/2(+1) nkx_plot=nkxy_plot/nx_kr
  call set_iymax(iymax,ix)
  do iy=0,iymax
   do ix3=0,nx_kr-1	!index of site inside kr cell
    do iy3=0,ny_kr-1
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)
    do icf=0,1
     do is=0,1
      do o1=0, l_orb(icf)-1
      ind=index_xyzcso_red(ix,iy,ix3,iy3,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
      indp= index_xyzcso(ix3,iy3,iz,icf,is,o1,nx_kr,ny_kr, nz_g0) 
!all z; x and y only in the scattering region in the middle
 do ix2=ixmin1,ixmax1	!index of x point, not of cell
  do iy2=iymin1,iymax1
   do iz2=0, nz_scf-1
    do icf2=0,1
     do is2=0,1
      do o2=0,l_orb(icf2)-1
      ind2=index_xyzcso(ix2-ixmin1,iy2-iymin1,iz2,icf2,is2,o2,nx, ny, nz_scf) 			!index 2 of G0_x (with V)
      indp2= index_xyzcso(mod(ix2,nx_kr),mod(iy2,ny_kr),iz2,icf2,is2,o2,nx_kr,ny_kr, nz_g0) 	!index_1 of G0_x (on the surface)


ind3=index_g0xtmp(ix-ix2/nx_kr,iy-iy2/ny_kr,nkx_plot,nky_plot)			!x and y of cell


!allocate(G0_x1(0:size_red-1,0:size_v_red-1))	!G_0(z,z',r,r') surface - scatt
!size_red=nxy_red*nz_plot*2*(l_orbital+l_orbital_f)*nx_kr*ny_kr
!size_v_red=2*nx*ny*nz_scf*(l_orbital+l_orbital_f)	!only the first nz_scf layers
!allocate(G0_x_tmp(0:nkx_plot*nky_plot-1,0:size_h_kr_red-1,0:size_h_kr_red-1)) !G_0(r-r',z,z')

G0_x1(ind, ind2)=G0_x_tmp(indp,indp2,ind3)		!surface - scatt

ind3=index_g0xtmp(-ix+ix2/nx_kr,-iy+iy2/ny_kr,nkx_plot,nky_plot)			!x and y of cell
G0_x2(ind2, ind)=G0_x_tmp(indp2,indp,ind3)		!scatt - surface



enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


!G0 for small system
!print *, "G0(r,r',z,z') small"

 do ix=0,nx-1		!index of site
  do iy=0,ny-1
   do iz=0,nz_scf-1	!z
    do icf=0,1
     do is=0,1
      do o1=0,l_orb(icf)-1
      ind=index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz_scf) !ix + iy*nx + iz*nx*nx + icf*nx*nx*nz_scf + is*2*nx*nx*nz_scf + o1*4*nx*nx*nz_scf
      indp=index_xyzcso(mod(ix+ixmin1,nx_kr),mod(iy+iymin1,ny_kr),iz,icf,is,o1,nx_kr,ny_kr,nz_g0)! iz + icf*nz_scf + is*2*nz_scf + o1*4*nz_scf
 do ix2=0,nx-1
  do iy2=0,ny-1
   do iz2=0,nz_scf-1	!z
    do icf2=0,1
     do is2=0,1
      do o2=0,l_orb(icf2)-1
      ind2=index_xyzcso(ix2,iy2,iz2,icf2,is2,o2,nx, ny, nz_scf)!ix2 + iy2*nx + iz2*nx*nx + icf2*nx*nx*nz_scf + is2*2*nx*nx*nz_scf + o2*4*nx*nx*nz_scf
      indp2=index_xyzcso(mod(ix2+ixmin1,nx_kr),mod(iy2+iymin1,ny_kr),iz2,icf2,is2,o2,nx_kr,ny_kr,nz_g0) !iz2 + icf2*nz_scf + is2*2*nz_scf + o2*4*nz_scf


ind3=index_g0xtmp((ix+ixmin1)/nx_kr-(ix2+ixmin1)/nx_kr,(iy+iymin1)/ny_kr-(iy2+iymin1)/ny_kr,nkx_plot,nky_plot)
!ind3=index_g0xtmp(ix2-ix,iy2-iy,nkx_plot,nx_red)
	
G0(ind, ind2)=G0_x_tmp(indp,indp2,ind3)		!scatt - scatt
!G0_x_b(ind, ind2)=G0_x_tmp_b(ind3,indp,indp2)


enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo




time_1=secnds(0.0)
time_g0rr=time_g0rr+time_1-time_0



endif

deallocate(G0_x_tmp)


time_0=secnds(0.0)

if (lgreen) then

print *, "G_0V"

 allocate(G0V(0:size_v_red-1,0:size_v_red-1))
 G0V=0
 call matmul1(size_v_red, G0, V_red, G0V)

!I-G0V
G0V=Id_red-G0V

print *, "(I-G0V)^-1"

 allocate(temp(0:size_v_red-1,0:size_v_red-1))
 temp=0

 call invertcmplx(size_v_red, G0V, temp, Id_red)		!AX=B (i omega - Ham) G1=Id

 deallocate(G0V)

print *, "T=V(I-G0V)^-1"

allocate(T_m(0:size_v_red-1,0:size_v_red-1))
T_m=0

 call matmul1(size_v_red, V_red, temp, T_m)

deallocate(temp)
!deallocate(G0)



if (nx_kr==nx .and. ny_kr==ny) then		!1 cell scattering

!WRITE(UNIT=ae2,FMT='(SPf7.4)') omega

!OPEN(unit=100,file=trim(label)//'/g/t_g_red_e'//ae2,status='unknown')

!do i=0,dim_hk-1

! write(100,"(10(a1,e10.3,a1,e10.3,a3))") ( "(",real(T_m(i,j)),",",aimag(T_m(i,j)),")  ", j=0,dim_hk-1 )

!enddo

! close(unit=100)

OPEN(unit=100,file=trim(label)//'/g/v_g_red',status='unknown')

do i=0,size_v_red-1

 write(100,"(16(a1,e8.1,a1,e8.1,a2))") ( "(",real(V_red(i,j)),",",aimag(V_red(i,j)),") ", j=0,size_v_red-1 )

enddo

 close(unit=100)

endif
!

if (lphaseshifts) then

print *, "Phase shifts"

!Phase shifts
!exp(2i phi_i)
!allocate(S_m(0:size_v_red-1,0:size_v_red-1))
!allocate(exphi(0:size_v_red-1))
allocate(exphi_t(0:size_v_red-1))
allocate(exphi_v(0:size_v_red-1))
exphi_t=0
exphi_v=0

!S_m=Id_red-2*pi*(0,1)*T_m

!call DIAGONALIZE_s(size_v_red, S_m, exphi )
call DIAGONALIZE_s(size_v_red, T_m, exphi_t )
call DIAGONALIZE_s(size_v_red, V_red,   exphi_v )


OPEN(unit=100,file=trim(label)//'/g/phase_shift',status='unknown', position='append')

do i=0, size_v_red-1
write(100,"(i5, f12.5, 40e15.6)"), i, omega, phase_of_t(exphi_t(i)),phase_of_t(exphi_v(i)), real(exphi_t(i))**2+aimag(exphi_t(i))**2, real(exphi_v(i))**2+aimag(exphi_v(i))**2
!real(exphi_t(i)), aimag(exphi_t(i)), real(exphi_v(i)), aimag(exphi_v(i))
!, real(exphi(i)), aimag(exphi(i)), (real(exphi(i)))**2+(aimag(exphi(i)))**2, phase_of_exp(exphi(i)) 
enddo

 close(unit=100)
!deallocate(s_m)
!deallocate(exphi)
deallocate(exphi_t)
deallocate(exphi_v)

endif



time_1=secnds(0.0)
time_gg0t=time_gg0t+time_1-time_0



print *, "g0T"
time_0=secnds(0.0)

allocate(G0T(0:size_red-1,0:size_v_red-1))
G0T=0

!call matmul(size_large, G0_x, T_m, G0T)
call matmul2(size_red,size_v_red,size_v_red, G0_x1, T_m, G0T)

deallocate(T_m)


print *, "g0Tg0"

!allocate(G_fin(0:size_large-1,0:size_large-1))
!allocate(G_fin(0:size_red-1,0:size_red-1))
allocate(dos_g(0:size_red-1))
allocate(dos_cf(0:nxy_red*nz_plot*nx_kr*ny_kr-1,0:dim_hk-1,0:dim_hk-1))
!G_fin=0
dos_g=0
dos_cf=0

!call matmul2(size_red,size_red,size_v_red, G0T, G0_x2, G_fin)
!call matmul1(size_large, G0T, G0_x, G_fin)

! do ix=ixmin1,ixmax1	!index of site
!  do iy=iymin1,iymax1
!   do iz=0, nz_scf-1
!    do icf=0,1
!     do is=0,1
!      do o1=0,l_orb(icf)-1
!      ind= index_xyzcso(ix-ixmin1,iy-iymin1,iz,icf,is,o1,nx, ny, nz_scf)

! do ix2=0,nx_red-1		!index of kr cell
! if (lc4v) iymax=ix2
! if (.not. lc4v) iymax=ny_red-1
!  do iy2=0,iymax
!   do ix4=0,nx_kr-1	!index of site inside kr cell
!    do iy4=0,ny_kr-1
!  do iz_ind=0,nz_plot-1	!z
!    iz2=z_values_plot(iz_ind)
!    do icf2=0,1
!     do is2=0,1
!      do o2=0,l_orb(icf2)-1
!      ind2=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,icf2,is2,o2,nz_plot,nx_red,nxy_red)

 call set_ixmax(ixmax)
 do ix=0,ixmax		!index of kr cell nx_red=nkx_plot/2(+1) nkx_plot=nkxy_plot/nx_kr
  call set_iymax(iymax,ix)
  do iy=0,iymax
   do ix3=0,nx_kr-1	!index of site inside kr cell
    do iy3=0,ny_kr-1
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)
    do icf=0,1
     do is=0,1
      do o1=0, l_orb(icf)-1
      ind=index_xyzcso_red(ix,iy,ix3,iy3,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
!      indp= index_xyzcso(ix3,iy3,iz,icf,is,o1,nx_kr,ny_kr, nz_scf) 
!all z; x and y only in the scattering region in the middle
 do ix2=ixmin1,ixmax1	!index of x point, not of cell
  do iy2=iymin1,iymax1
   do iz2=0, nz_scf-1
    do icf2=0,1
     do is2=0,1
      do o2=0,l_orb(icf2)-1
      ind2=index_xyzcso(ix2-ixmin1,iy2-iymin1,iz2,icf2,is2,o2,nx, ny, nz_scf) 			!index 2 of G0_x (with T) and of G0T
!      indp2= index_xyzcso(mod(ix2,nx_kr),mod(iy2,ny_kr),iz2,icf2,is2,o2,nx_kr,ny_kr, nz_scf) 	!

dos_g(ind)=dos_g(ind)+G0T(ind,ind2)*G0_x2(ind2,ind)

    do icf3=0,1
     do is3=0,1
      do o3=0,l_orb(icf3)-1
      ind3=index_xyzcso_red(ix,iy,ix3,iy3,iz_ind,icf3,is3,o3,nz_plot,nx_red,nxy_red)	!same position as ind, but different orbitals
      indp=index_xyz_red(ix,iy,ix3,iy3,iz_ind,nz_plot, nx_red, nxy_red)	!position of atom
      ind2p=index_cso(icf,is,o1) !icf2+2*is2+4*o2
      ind3p=index_cso(icf3,is3,o3) !icf3+2*is3+4*o3

dos_cf(indp,ind2p,ind3p)= dos_cf(indp,ind2p,ind3p)+G0T(ind,ind2)*G0_x2(ind2,ind3)

     enddo
    enddo
    enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


!print *, "dos_g"
!ix2=nx/2
!iy2=ix2
!iz_ind=0
!indp=ix2 + iy2*nx_red - (iy2*(iy2+1))/2 + iz_ind*nxy_red

!do i=0,dim_hk-1
!write (*, "(10(1x,f8.4,f8.4))") (dos_cf(indp, i,j), j=0,dim_hk-1)
!enddo
!print *, ""

!if (l_orbital==2) then
!write (*, "(8(1x,a16))") "c1-", "f1-", "c1+", "f1+", "c2-", "f2-", "c2+", "f2+"
!do i=0,7
!write (*, "(8(1x,f8.4,f8.4))") (dos_cf(0, i,j), j=0,7)
!enddo
!endif
!print *, ""
!write (*, "(8(1x,f8.4,f8.4))") (dos_g(j), j=0,7)





deallocate(G0T)

!G_fin=G_fin+G0_x
!dos_g=dos_g+dos_0

endif




if (lborn) then
print *, "g0V"

allocate(G0T(0:size_red-1,0:size_v_red-1))
G0T=0

!call matmul1(size_large, G0_x, T_m, G0T)
call matmul2(size_red,size_v_red,size_v_red, G0_x1, V_red, G0T)			!V instead of T


print *, "g0Vg0"

!allocate(G_fin(0:size_large-1,0:size_large-1))
!allocate(G_born(0:size_red-1,0:size_red-1))
!G_born=0
allocate(dos_born(0:size_red-1))
allocate(dos_cf_born(0:nxy_red*nz_plot*nx_kr*ny_kr-1,0:dim_hk-1,0:dim_hk-1))
dos_born=0
dos_cf_born=0

!call matmul2(size_red,size_red,size_v_red, G0T, G0_x2, G_born)

 call set_ixmax(ixmax)
 do ix=0,ixmax	!index of kr cell nx_red=nkx_plot/2(+1) nkx_plot=nkxy_plot/nx_kr
  call set_iymax(iymax,ix)
  do iy=0,iymax
   do ix3=0,nx_kr-1	!index of site inside kr cell
    do iy3=0,ny_kr-1
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)
    do icf=0,1
     do is=0,1
      do o1=0, l_orb(icf)-1
      ind=index_xyzcso_red(ix,iy,ix3,iy3,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
!      indp= index_xyzcso(ix3,iy3,iz,icf,is,o1,nx_kr,ny_kr, nz_scf) 
!all z; x and y only in the scattering region in the middle
 do ix2=ixmin1,ixmax1	!index of x point, not of cell
  do iy2=iymin1,iymax1
   do iz2=0, nz_scf-1
    do icf2=0,1
     do is2=0,1
      do o2=0,l_orb(icf2)-1
      ind2=index_xyzcso(ix2-ixmin1,iy2-iymin1,iz2,icf2,is2,o2,nx, ny, nz_scf) 			!index 2 of G0_x (with T) and of G0T
!      indp2= index_xyzcso(mod(ix2,nx_kr),mod(iy2,ny_kr),iz2,icf2,is2,o2,nx_kr,ny_kr, nz_scf) 	!

dos_born(ind)=dos_born(ind)+G0T(ind,ind2)*G0_x2(ind2,ind)

    do icf3=0,1
     do is3=0,1
      do o3=0,l_orb(icf3)-1
      ind3=index_xyzcso_red(ix,iy,ix3,iy3,iz_ind,icf3,is3,o3,nz_plot,nx_red,nxy_red)	!same position as ind, but different orbitals
      indp=index_xyz_red(ix,iy,ix3,iy3,iz_ind,nz_plot, nx_red, nxy_red)	!position of atom
      ind2p=index_cso(icf,is,o1) !icf2+2*is2+4*o2
      ind3p=index_cso(icf3,is3,o3) !icf3+2*is3+4*o3

dos_cf_born(indp,ind2p,ind3p)= dos_cf_born(indp,ind2p,ind3p)+G0T(ind,ind2)*G0_x2(ind2,ind3)

     enddo
    enddo
    enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo



!print *, "dos_g_born"
!ix2=nx/2
!iy2=ix2
!iz_ind=0
!indp=ix2 + iy2*nx_red - (iy2*(iy2+1))/2 + iz_ind*nxy_red

!do i=0,dim_hk-1
!write (*, "(10(1x,f8.4,f8.4))") (dos_cf_born(indp, i,j), j=0,dim_hk-1)
!enddo
!print *, ""



deallocate(G0T)

!G_born=G_born+G0_x
!dos_born=dos_born+dos_0

endif


deallocate(G0)
!deallocate(G0_x)
deallocate(G0_x1)
deallocate(G0_x2)

time_1=secnds(0.0)
time_g0tg0=time_g0tg0+time_1-time_0



!!!!!!!!!!!!!!!!!!!Final DOS
print *, "DOS"

time_0=secnds(0.0)

entry_dos=18
allocate(dos_fin(0:nkx_plot*nx_kr-1,0:nky_plot*ny_kr-1,0:nz_plot-1,entry_dos))
dos_fin=0
if(all_g) allocate(g_fin(0:nkx_plot*nx_kr-1,0:nky_plot*ny_kr-1,0:nz_plot-1,2*dim_hk**2))
if(all_g) g_fin=0
!1: rho_d1
!2: rho_d2
!3: rho_f1
!4: rho_f2
!5: rho_f3
!6: -1/pi Im Gd1f1
!7: -1/pi Im Gd1f3
!8: -1/pi Im Gf1f3
!9: -1/pi Im Gd2f2
!10-18 same for G0 
!19:27 same for G_born


!dos_fin indexes
do ix=0, nkx_plot-1	!index of kr cell
 do iy=0, nky_plot-1	
  do ix3=0, nx_kr-1	!index of atom inside kr cell
   do iy3=0, ny_kr-1
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)

ix2=ix	!index of cell
iy2=iy
ix4=ix3	!index of atom inside the cell
iy4=iy3

!G indexes
if(lc4v) then		!nx_kr=1, ny_kr=1
if(ix < nx_red) ix2=ix
if(ix > nx_red-1) ix2=nkx_plot-1-ix
if(iy < nx_red) iy2=iy
if(iy > nx_red-1) iy2=nkx_plot-1-iy
!if    (ix < nx_red .and. iy <=ix) then
!ix2=ix
!iy2=iy
!elseif(ix < nx_red .and. iy > ix) then
!ix2=iy
!iy2=ix
!elseif(ix > nx_red-1 .and. iy > ix) then
if (iy2>ix2) then
ix_temp=ix2
ix2=iy2
iy2=ix_temp
endif

else if(lcs) then
if(iy > ny_red-1) iy2=nky_plot-1-iy 	!mirror symmetry

else if(lnos) then
!do nothing
endif



ix5=ix3+ix*nx_kr
iy5=iy3+iy*ny_kr
!write(*, "(20i5)"), ix,iy,ix2,iy2,ix3,iy3,ix4,iy4,ix5,iy5

!1:(G-G0)_c
!2:(G-G0)_f
!3:G0_c
!4:G0_f
!5:(G-G0)_c born
!6:(G-G0)_f born
!7-12, (real parts?, magnetic response?) tGt

!rho_c, summed over spin and orbital
!      ind=index_xyzcso_red(ix ,iy ,ix2,iy2,iz_ind,icf,is,o1,nz_plot, nx_red, nxy_red)
      ind1=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,0,0,0,nz_plot,nx_red,nxy_red)	!c1d
      ind2=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,0,1,0,nz_plot,nx_red,nxy_red)	!c1u
      ind3=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,0,0,1,nz_plot,nx_red,nxy_red)	!c2d
      ind4=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,0,1,1,nz_plot,nx_red,nxy_red) 	!c2u

if (lgreen .and. l_orbital==1) dos_fin(ix5,iy5,iz,1)=-aimag(dos_g(ind1)+dos_g(ind2))/pi		!G-G0+G0
if (lgreen .and. l_orbital==2) dos_fin(ix5,iy5,iz,1)=-aimag(dos_g(ind1)+dos_g(ind2)+dos_g(ind3)+dos_g(ind4))/pi		!G-G0+G0

if (l_orbital==1)             dos_fin(ix5,iy5,iz,3)=-aimag(dos_0(ind1)+dos_0(ind2))/pi			!G0
if (l_orbital==2)             dos_fin(ix5,iy5,iz,3)=-aimag(dos_0(ind1)+dos_0(ind2)+dos_0(ind3)+dos_0(ind4))/pi			!G0

if (lborn .and. l_orbital==1)  dos_fin(ix5,iy5,iz,5)= -aimag(dos_born(ind1)+dos_born(ind2))/pi	!G_born-G0+G0
if (lborn .and. l_orbital==2)  dos_fin(ix5,iy5,iz,5)= -aimag(dos_born(ind1)+dos_born(ind2)+dos_born(ind3)+dos_born(ind4))/pi	!G_born-G0+G0

!rho_f, summed over spin and orbital
      ind1=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,0,0,nz_plot,nx_red,nxy_red)	!f1d 
      ind2=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,1,0,nz_plot,nx_red,nxy_red) 	!f1u
      ind3=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,0,1,nz_plot,nx_red,nxy_red) 	!f2d
      ind4=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,1,1,nz_plot,nx_red,nxy_red) 	!f2u
      ind5=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,0,2,nz_plot,nx_red,nxy_red) !gamma7!f3d
      ind6=index_xyzcso_red(ix2,iy2,ix4,iy4,iz_ind,1,1,2,nz_plot,nx_red,nxy_red) !gamma7!f3u

if (lgreen .and. l_orbital_f==1) dos_fin(ix5,iy5,iz,2)=-aimag(dos_g(ind1)+dos_g(ind2))/pi	!G
if (lgreen .and. l_orbital_f==2) dos_fin(ix5,iy5,iz,2)=-aimag(dos_g(ind1)+dos_g(ind2)+dos_g(ind3)+dos_g(ind4))/pi	!G
if (lgreen .and. l_orbital_f==3) dos_fin(ix5,iy5,iz,2)=-aimag(dos_g(ind1)+dos_g(ind2)+dos_g(ind3)+dos_g(ind4)+dos_g(ind5)+dos_g(ind6))/pi	!G

if (l_orbital_f==1)            dos_fin(ix5,iy5,iz,4)= -aimag(dos_0(ind1)+dos_0(ind2))/pi	!G0
if (l_orbital_f==2)            dos_fin(ix5,iy5,iz,4)= -aimag(dos_0(ind1)+dos_0(ind2)+dos_0(ind3)+dos_0(ind4))/pi	!G0
if (l_orbital_f==3)            dos_fin(ix5,iy5,iz,4)= -aimag(dos_0(ind1)+dos_0(ind2)+dos_0(ind3)+dos_0(ind4)+dos_0(ind5)+dos_0(ind6))/pi	!G0

if (lborn .and. l_orbital_f==1)  dos_fin(ix5,iy5,iz,6)= -aimag(dos_born(ind1)+dos_born(ind2))/pi !G_born
if (lborn .and. l_orbital_f==2)  dos_fin(ix5,iy5,iz,6)= -aimag(dos_born(ind1)+dos_born(ind2)+dos_born(ind3)+dos_born(ind4))/pi !G_born
if (lborn .and. l_orbital_f==3)  dos_fin(ix5,iy5,iz,6)= -aimag(dos_born(ind1)+dos_born(ind2)+dos_born(ind3)+dos_born(ind4)+dos_born(ind5)+dos_born(ind6))/pi !G_born

![tGt]
!ind5p=ix2 + iy2*nx_red - (iy2*(iy2+1))/2 + iz_ind*nxy_red
 ind5p=index_xyz_red(ix2,iy2,ix4,iy4,iz_ind,nz_plot,nx_red,nxy_red)
!ind5p=index_xyz_red(ix ,iy ,ix2,iy2,iz_ind,nz_plot,nx_red,nxy_red)	!position of atom

do i=0,dim_hk-1
do j=0,dim_hk-1
if (lgreen) 	dos_fin(ix5,iy5,iz,7)=dos_fin(ix5,iy5,iz,7)-aimag(t_hop(i)*dos_cf(ind5p,i,j)*t_hop(j))/pi
		dos_fin(ix5,iy5,iz,8)=dos_fin(ix5,iy5,iz,8)-aimag(t_hop(i)*dos_cf_0(ind5p,i,j)*t_hop(j))/pi
if (lborn) 	dos_fin(ix5,iy5,iz,9)=dos_fin(ix5,iy5,iz,9)-aimag(t_hop(i)*dos_cf_born(ind5p,i,j)*t_hop(j))/pi

!only for reconstructed surface:
!first layer: signal from existing row
!second layer: signal under non-existing row rescaled by alfa_nz
!z independent, for simplicity put only in iz=0
if (nx_kr==2 .and. lrec .and. iz==0 .and. nz_plot>1) then	!iz==0 so I avoid double counting
 ind6p=index_xyz_red(ix2,iy2,ix4,iy4,0,nz_plot,nx_red,nxy_red)
 ind7p=index_xyz_red(ix2,iy2,ix4,iy4,1,nz_plot,nx_red,nxy_red)
if (lgreen .and. mod(ix5,2)==0) 	dos_fin(ix5,iy5,0,10)=dos_fin(ix5,iy5,0,10)-aimag(t_hop(i)*dos_cf(ind6p,i,j)*t_hop(j))/pi
if (lgreen .and. mod(ix5,2)==1) 	dos_fin(ix5,iy5,0,10)=dos_fin(ix5,iy5,0,10)-aimag(t_hop(i)*dos_cf(ind7p,i,j)*t_hop(j))/pi*alfa_nz
if (mod(ix5,2)==0)			dos_fin(ix5,iy5,0,11)=dos_fin(ix5,iy5,0,11)-aimag(t_hop(i)*dos_cf_0(ind6p,i,j)*t_hop(j))/pi
if (mod(ix5,2)==1)			dos_fin(ix5,iy5,0,11)=dos_fin(ix5,iy5,0,11)-aimag(t_hop(i)*dos_cf_0(ind7p,i,j)*t_hop(j))/pi*alfa_nz
if (lborn .and. mod(ix5,2)==0) 		dos_fin(ix5,iy5,0,12)=dos_fin(ix5,iy5,0,12)-aimag(t_hop(i)*dos_cf_born(ind6p,i,j)*t_hop(j))/pi
if (lborn .and. mod(ix5,2)==1) 		dos_fin(ix5,iy5,0,12)=dos_fin(ix5,iy5,0,12)-aimag(t_hop(i)*dos_cf_born(ind7p,i,j)*t_hop(j))/pi*alfa_nz

!here for the DOS
do icf=0,1
do is=0,1
do o1=0, l_orb(icf)-1
      ind1=index_xyzcso_red(ix2,iy2,ix4,iy4,0,icf,is,o1,nz_plot,nx_red,nxy_red)	
      ind2=index_xyzcso_red(ix2,iy2,ix4,iy4,1,icf,is,o1,nz_plot,nx_red,nxy_red)	
if (lgreen .and. mod(ix5,2)==0) 	dos_fin(ix5,iy5,0,16)=dos_fin(ix5,iy5,0,16)-aimag(dos_g(ind1))/pi
if (lgreen .and. mod(ix5,2)==1) 	dos_fin(ix5,iy5,0,16)=dos_fin(ix5,iy5,0,16)-aimag(dos_g(ind2))/pi*alfa_nz
if (mod(ix5,2)==0)			dos_fin(ix5,iy5,0,17)=dos_fin(ix5,iy5,0,17)-aimag(dos_0(ind1))/pi
if (mod(ix5,2)==1)			dos_fin(ix5,iy5,0,17)=dos_fin(ix5,iy5,0,17)-aimag(dos_0(ind2))/pi*alfa_nz
if (lborn .and. mod(ix5,2)==0) 		dos_fin(ix5,iy5,0,18)=dos_fin(ix5,iy5,0,18)-aimag(dos_born(ind1))/pi
if (lborn .and. mod(ix5,2)==1) 		dos_fin(ix5,iy5,0,18)=dos_fin(ix5,iy5,0,18)-aimag(dos_born(ind2))/pi*alfa_nz
enddo
enddo
enddo
endif
enddo !i
enddo !j

!all_g
if(all_g) then
!print *, "Doing all_g"
!ind5p=ix2 + iy2*nx_red - (iy2*(iy2+1))/2 + iz_ind*nxy_red
ind5p=index_xyz_red(ix2,iy2,ix4,iy4,iz_ind,nz_plot,nx_red,nxy_red)

do i=0,dim_hk-1
do j=0,dim_hk-1
 if(lgreen) g_fin(ix,iy,iz,i+j*dim_hk+1)=		-aimag(dos_cf(ind5p,i,j))/pi	
 if(lborn)  g_fin(ix,iy,iz,i+j*dim_hk+1+dim_hk**2)=	-aimag(dos_cf_born(ind5p,i,j))/pi	
enddo
enddo


endif




enddo
enddo
  enddo
 enddo
enddo


!here I take the average on neighbouring columns
if (nx_kr==2 .and. lrec) then	!iz==0 so I avoid double counting
iz=0
do ix=0, nkx_plot-1	!index of kr cell
 do iy=0, nky_plot-1	

ix3=0
iy3=0
ix5=ix3+ix*nx_kr
iy5=iy3+iy*ny_kr

dos_fin(ix5,iy5,0,13)=dos_fin(ix5,iy5,0,7)	!on row=same as 9
dos_fin(ix5,iy5,0,14)=dos_fin(ix5,iy5,0,8)	!on row=same as 10
dos_fin(ix5,iy5,0,15)=dos_fin(ix5,iy5,0,9)	!on row=same as 11

enddo
enddo


do ix=0, nkx_plot-1	!index of kr cell
 do iy=0, nky_plot-1	
iz=0

ix3=1
iy3=0
ix5=ix3+ix*nx_kr
iy5=iy3+iy*ny_kr

dos_fin(ix5,iy5,0,13)=(dos_fin(mod(ix5-1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,7)+dos_fin(mod(ix5+1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,7))*stm_layer_ratio/2
dos_fin(ix5,iy5,0,14)=(dos_fin(mod(ix5-1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,8)+dos_fin(mod(ix5+1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,8))*stm_layer_ratio/2
dos_fin(ix5,iy5,0,15)=(dos_fin(mod(ix5-1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,9)+dos_fin(mod(ix5+1+nkx_plot*nx_kr,nkx_plot*nx_kr),iy5,0,9))*stm_layer_ratio/2

enddo
enddo
endif

!if (nz_plot>=2)  dos_fin(:,:,1,:)=dos_fin(:,:,1,:)*stm_layer_ratio	!rescales second layer signal by stm_layer_ratio needed in nztot
!if (nz_plot==3)  dos_fin(:,:,2,:)=dos_fin(:,:,2,:)*stm_layer_ratio**2	!rescales third layer signal by stm_layer_ratio**2



time_1=secnds(0.0)
time_dos=time_dos+time_1-time_0


!Writes to file
time_0_w=secnds(0.0)

if (lqpi) then

!WRITE(UNIT=ae2,FMT='(SPf7.4)') omega

!stm spatial map, fixed energy
!1,2:x,y  3,4:DeltaG, 5,6:G0. 7,8:DeltaG_born 9,10,11:tGt

do iz_ind=0,nz_plot-1	!z
 iz=z_values_plot(iz_ind)
 WRITE(UNIT=az,FMT='(I2.2)') iz
 OPEN(unit=100,file=trim(label)//'/g/dos_nz'//az//'_e'//ae2,status='unknown')

 do ix=0, nkx_plot*nx_kr-1
  do iy=0, nky_plot*ny_kr-1

!write(100,"(2i5,30e14.5)"),  ix,iy,(dos_fin(ix,iy,iz,i), i=1,entry_dos )
write(100,"(2f10.5,30e14.5)"),  x_ss(ix,iy),y_ss(ix,iy),(dos_fin(ix,iy,iz,i), i=1,entry_dos )
!				dos_fin(ix,iy,iz,1),dos_fin(ix,iy,iz,2),&	!DeltaG_c, DeltaG_f
!				 dos_fin(ix,iy,iz,3), dos_fin(ix,iy,iz,4),&	!G0_c,G0_f
!				  dos_fin(ix,iy,iz,5), dos_fin(ix,iy,iz,6), &	!G_c, G_f born
!				  dos_fin(ix,iy,iz,7), dos_fin(ix,iy,iz,8), dos_fin(ix,iy,iz,9),&!, tGt, tG0t, tGbornt
!				  dos_fin(ix,iy,iz,10), dos_fin(ix,iy,iz,11), dos_fin(ix,iy,iz,12)!, tGt, tG0t, tGbornt alternatively 1st and secondlayer


  enddo	!ix
 enddo	!iy
 write(100,"(a5,f15.6)") "#T_g=", T_g
 write(100,"(a5,10f15.6)") "#t=", (t_hop(i), i=0, dim_hk-1)
 close(unit=100)

enddo	!iz


!summed over z
! OPEN(unit=100,file=trim(label)//'/g/dos_nztot_e'//ae2,status='unknown')
! do ix=0, nkx_plot*nx_kr-1
!  do iy=0, nky_plot*ny_kr-1

!write(100,"(2i5,30e14.5)"), ix,iy,(sum(dos_fin(ix,iy,:,j)),j=1,entry_dos)

!  enddo
! enddo
! write(100,"(a5,f15.6)") "#T_g=", T_g
! write(100,"(a5,10f15.6)") "#t=", (t_hop(i), i=0, dim_hk-1)
! close(unit=100)


 
 endif	!lqpi
!write(100,*), ""

!print *, "DOS print"




!stm energy scan, fixed position, only close  to the center
!1,2:x,y  3,4:Gtot, 5,6:G0. 7,8:G_born
if (mod(nkx_plot*nx_kr,2)==1) ixcenter=(nkx_plot*nx_kr)/2
if (mod(nky_plot*ny_kr,2)==1) iycenter=(nky_plot*ny_kr)/2
if (mod(nkx_plot*nx_kr,2)==0) ixcenter=(nkx_plot*nx_kr)/2-1
if (mod(nky_plot*ny_kr,2)==0) iycenter=(nky_plot*ny_kr)/2-1


ixplotmin=max(ixcenter-2,0)
iyplotmin=max(iycenter-2,0)
if(lc4v .or. lc2v) ixplotmax=ixcenter
if(lcs  .or. lnos) ixplotmax=min(ixcenter+2,nkx_plot*nx_kr-1) 

!print *, "ixcenter,iycenter:", ixcenter, iycenter
!print *, "ixplotmin,ixplotmax:", ixplotmin, ixplotmax

 do ix=ixplotmin,ixplotmax	!x
if(lc4v) iyplotmax=ix
if(lc2v .or. lcs) iyplotmax=iycenter
if(lnos) iyplotmax=min(iycenter+2,nky_plot*ny_kr-1) 
!print *, "ix,iyplotmin,iyplotmax", ix,iyplotmin, iyplotmax
  do iy=iyplotmin,iyplotmax	!y
   do iz_ind=0,nz_plot-1	!z
    iz=z_values_plot(iz_ind)

    WRITE(UNIT=ax4,FMT='(I4.4)') ix
    WRITE(UNIT=ay4,FMT='(I4.4)') iy
    WRITE(UNIT=az,FMT='(I2.2)') iz
 
    ind=ix +iy*nx_red +iz*nx_red*nx_red
    unitfile=200+ind
    OPEN(unit=unitfile,file=trim(label)//'/g/dos_nx'//ax4//'_ny'//ay4//'_nz'//az,status='unknown', position='append')

    write(unitfile,"(2f10.5,30e14.5)"), omega,omega, (dos_fin(ix,iy,iz,j), j=1,entry_dos)! &! (dos_fin(ix,iy,iz,i)+dos_fin(ix,iy,iz,i+10), i=1,10),&	!G+G0
    						!(dos_fin(ix,iy,iz,i+10), i=1,10),&			!G0
    						!(dos_fin(ix,iy,iz,i+20)+dos_fin(ix,iy,iz,i+10), i=1,10)	!Gborn+G0
!    				dos_fin(ix,iy,iz,1)+dos_fin(ix,iy,iz,3), dos_fin(ix,iy,iz,2)+dos_fin(ix,iy,iz,4),&	!G_c, G_f
!				 dos_fin(ix,iy,iz,3), dos_fin(ix,iy,iz,4),&								!G0_c, G0_f
!				  dos_fin(ix,iy,iz,5)+dos_fin(ix,iy,iz,3), dos_fin(ix,iy,iz,6)+dos_fin(ix,iy,iz,4), &			!G_c born, G_f born
!				  dos_fin(ix,iy,iz,7)+dos_fin(ix,iy,iz,8), dos_fin(ix,iy,iz,8),dos_fin(ix,iy,iz,9)+dos_fin(ix,iy,iz,8),&	!tGt tG0t tGbornt
!				  dos_fin(ix,iy,iz,10)+dos_fin(ix,iy,iz,11), dos_fin(ix,iy,iz,11),dos_fin(ix,iy,iz,12)+dos_fin(ix,iy,iz,11)	!tGt tG0t tGbornt
				

    close(unit=unitfile)

if(all_g) then
    OPEN(unit=unitfile,file=trim(label)//'/g/g_nx'//ax4//'_ny'//ay4//'_nz'//az,status='unknown', position='append')

!ind5p=ix + iy*nx_red - (iy*(iy+1))/2 + iz_ind*nxy_red
ind5p=index_xyz_red(ix/nx_kr,iy/ny_kr,mod(ix2,nx_kr),mod(iy2,ny_kr),iz_ind,nz_plot, nx_red, nxy_red)

if (lgreen .and. lborn)    write(unitfile,"(2f10.5,300e13.5)"), omega,omega, &
					 (-aimag(dos_cf(ind5p,mod(i,dim_hk),i/dim_hk)+dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1), &
					 (-aimag(dos_cf_born(ind5p,mod(i,dim_hk),i/dim_hk)+dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1),&
					 (-aimag(dos_cf_0(0,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1)					 
if (lgreen .and. .not. lborn)    write(unitfile,"(2f10.5,300e13.5)"), omega,omega, &
					 (-aimag(dos_cf(ind5p,mod(i,dim_hk),i/dim_hk)+dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1), &
					 (0., i=0, dim_hk**2-1),&
					 (-aimag(dos_cf_0(0,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1)					 
if (.not. lgreen .and. lborn)    write(unitfile,"(2f10.5,300e13.5)"), omega,omega, &
					 (0., i=0, dim_hk**2-1), &
					 (-aimag(dos_cf_born(ind5p,mod(i,dim_hk),i/dim_hk)+dos_cf_0(ind5p,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1),&
					 (-aimag(dos_cf_0(0,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1)					 
!if (.not. lgreen .and. .not. lborn)    write(unitfile,"(2f10.5,300e13.5)"), omega,omega, &
!					 (0., i=0, dim_hk**2-1), &
!					 (0., i=0, dim_hk**2-1), &
!					 (-aimag(dos_cf_0(0,mod(i,dim_hk),i/dim_hk))/pi, i=0, dim_hk**2-1)					 

!    write(*,"(2f10.5,200i4)"), omega,omega,  (mod(i,dim_hk),i/dim_hk, i=0, dim_hk**2-1)


    close(unit=unitfile)

endif

  enddo
 enddo
enddo


do ix=ixplotmin,ixplotmax	!x
 if(lc4v) iyplotmax=ix
 if(lc2v .or. lcs) iyplotmax=iycenter
 if(lnos) iyplotmax=min(iycenter+2,nky_plot*ny_kr-1) 
  do iy=iyplotmin,iyplotmax	!y

    WRITE(UNIT=ax4,FMT='(I4.4)') ix
    WRITE(UNIT=ay4,FMT='(I4.4)') iy
 
    ind=ix +iy*nx_red
    unitfile=200+ind
    OPEN(unit=unitfile,file=trim(label)//'/g/dos_nx'//ax4//'_ny'//ay4//'_nztot',status='unknown', position='append')

    write(unitfile,"(2f10.5,30e14.5)"), omega,omega, (sum(dos_fin(ix,iy,:,j)), j=1,entry_dos)
    !sum(dos_fin(ix,iy,:,1))+sum(dos_fin(ix,iy,:,3)), sum(dos_fin(ix,iy,:,2))+sum(dos_fin(ix,iy,:,4)),&	!G_c, G_f
	!			 sum(dos_fin(ix,iy,:,3)), sum(dos_fin(ix,iy,:,4)),&								!G0_c, G0_f
	!			  sum(dos_fin(ix,iy,:,5))+sum(dos_fin(ix,iy,:,3)), sum(dos_fin(ix,iy,:,6))+sum(dos_fin(ix,iy,:,4)), &			!G_c born, G_f born
	!			  sum(dos_fin(ix,iy,:,7))+sum(dos_fin(ix,iy,:,8)), sum(dos_fin(ix,iy,:,8)),sum(dos_fin(ix,iy,:,9))+sum(dos_fin(ix,iy,:,8)),&	!tGt tG0t tGbornt
	!			  sum(dos_fin(ix,iy,:,10))+sum(dos_fin(ix,iy,:,11)), sum(dos_fin(ix,iy,:,11)),sum(dos_fin(ix,iy,:,12))+sum(dos_fin(ix,iy,:,11))	!tGt tG0t tGbornt
enddo
enddo
				

    close(unit=unitfile)




time_1_w=secnds(0.0)
time_write=time_Write+time_1_w-time_0_w


deallocate(dos_0)
deallocate(dos_cf_0)
deallocate(dos_cf_r)
!deallocate(G0_x)
!deallocate(G0_x_b)
if (lgreen) deallocate(dos_g)
if (lgreen) deallocate(dos_cf)
!if (lborn) deallocate(G_born)
if (lborn) deallocate(dos_born)
if (lborn) deallocate(dos_cf_born)
!deallocate(G)
!deallocate(G_b)


!!!!!!!!!!!!Fourier transform of dos/g == QPI

!I remove G0, so only the signal is transformed (no peak at q=0)

time_fft_0=secnds(0.0)

if (lqpi .and. (lborn .or. lgreen)) then
print *, "Fourier"


l_dfti(1)=nkx_plot*nx_kr
l_dfti(2)=nky_plot*ny_kr


xcenter=(nkx_plot*nx_kr-1)/2	!point around which center the fourier transform (important to have it real)
ycenter=(nky_plot*ny_kr-1)/2
print *, "center:x,y=", xcenter, ycenter
!if (mod(nkx_plot*nx_kr,2)==0) exp_corrx_kr=0.5d0
!if (mod(nkx_plot*nx_kr,2)==1) exp_corrx_kr=0.5d0*dble(nkx_plot*nx_kr-1)/dble(nkx_plot*nx_kr)
!if (mod(nky_plot*ny_kr,2)==0) exp_corry_kr=0.5d0
!if (mod(nky_plot*ny_kr,2)==1) exp_corry_kr=0.5d0*dble(nky_plot*ny_kr-1)/dble(nky_plot*ny_kr)

if(all_g_fourier) allocate(g_fourier(0:nkr_plot*nx_kr*ny_kr-1,0:nz_plot-1,2*dim_hk**2))
if(all_g_fourier) g_fourier=0
allocate(dos_fourier(0:nkr_plot*nx_kr*ny_kr-1,0:nz_plot-1,entry_dos))
dos_fourier=0


allocate(dos_r_fft(0:nkx_plot*nky_plot*nx_kr*ny_kr-1))


do iz_ind=0,nz_plot-1	
 iz=z_values_plot(iz_ind)
 do i=1,entry_dos
 
  dos_r_fft=0


  do ix=0,nkx_plot*nx_kr-1
   do iy=0,nky_plot*ny_kr-1

    ind=ix + iy*nkx_plot*nx_kr
    dos_r_fft(ind)=dos_fin(ix,iy,iz_ind,i)*&
    			exp(-2*pi*(0,1)*((ix-xcenter)*exp_corrx_kr+(iy-iycenter)*exp_corry_kr))!corr fct to center around (n-1)/2n in x space
!    			exp(2*pi*(0,1)*(exp_corrx_kr*ix+exp_corry_kr*iy))!2nd option bugged
!    			exp(2*pi*(0,1)*(ix+iy-(nkx_plot-1))*(nkx_plot-1)/2/nkx_plot)!corr fct to center around (n-1)/2n in x space
   enddo
  enddo

l_dfti(1)=nkx_plot*nx_kr
l_dfti(2)=nky_plot*ny_kr

!  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

!  status = DftiCommitDescriptor(hand)
! if (0 /= status) print *, "error"

!   status = DftiComputeForward(hand, dos_r_fft)		!-
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"

  do ix=0,nkx_plot*nx_kr-1
   do iy=0,nky_plot*ny_kr-1
    ind=ix + iy*nkx_plot*nx_kr
   dos_fourier(ind,iz_ind,i)=dos_r_fft(ind)*&
    			exp(2*pi*(0,1)*((dble(ix)/(nkx_plot*nx_kr)*ixcenter+(dble(iy)/(nky_plot*ny_kr)*iycenter))))
!    			exp(2*pi*(0,1)*((ix*exp_corrx_kr+iy*exp_corry)))	!ok only for odd
!    			exp(2*pi*(0,1)*((exp_corrx_kr-dble(ix)/(nkx_plot*nx_kr))*xcenter+(exp_corry-dble(iy)/(nky_plot*ny_kr))*ycenter))		!2nd option bugged
   		  	     !exp(2*(0,1)*pi*(ix+iy)*(nkx_plot-1)/(2*nkx_plot)) !corr fct to center around 0 for k vectors
   enddo
  enddo

 enddo

!!!all g
if(all_g_fourier) then
!print *, "Doing all_g fourier"

 do i=1, 2*dim_hk**2
 
  dos_r_fft=0
  do ix=0,nkx_plot*nx_kr-1
   do iy=0,nky_plot*ny_kr-1
    ind=ix + iy*nkx_plot*nx_kr
    dos_r_fft(ind)=g_fin(ix,iy,iz_ind,i)*&
    			exp(2*pi*(0,1)*(ix+iy-(nkx_plot-1))*(nkx_plot-1)/2/nkx_plot)!corr fct to center around (n-1)/2n in x space
   enddo
  enddo

!  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,l_dfti)
! if (0 /= status) print *, "error"

!  status = DftiCommitDescriptor(hand)
! if (0 /= status) print *, "error"

!  status = DftiComputeForward(hand, dos_r_fft)	!-
! if (0 /= status) print *, "error"

!  status = DftiFreeDescriptor(hand)
! if (0 /= status) print *, "error"

  do ix=0,nkx_plot*nx_kr-1
   do iy=0,nkx_plot*ny_kr-1
    ind=ix + iy*nkx_plot*nx_kr
   g_fourier(ind,iz_ind,i)=dos_r_fft(ind)*&
   		  	     exp(2*(0,1)*pi*(ix+iy)*(nkx_plot-1)/(2*nkx_plot)) !corr fct to center around 0 for k vectors
   enddo
  enddo

 enddo

endif	!all_g



enddo

deallocate(dos_r_fft)

time_1=secnds(0.0)
time_fourier=time_fourier+time_1-time_0

time_fft_1=secnds(0.0)
time_fft=time_fft+time_fft_1-time_fft_0

lfftreal=.true.
do iz_ind=0,nz_plot-1	
 iz=z_values_plot(iz_ind)
 do i=1,entry_dos
 
  do ix=0,nkx_plot*nx_kr-1
   do iy=0,nky_plot*ny_kr-1
    ind=ix + iy*nkx_plot*nx_kr
if   (abs(aimag(dos_fourier(ind,iz_ind,i)))>1e-7) lfftreal=.false.
enddo
enddo
enddo
enddo

print *, "is FFT real?", lfftreal

!Writes to file fourier

!momentum map, fixed energy
!1,2: kx,ky 3,4:G-G0, 5,6: G0, 7,8: G_born -G_0, 9,10,11:tGt

time_0_w=secnds(0.0)

do iz_ind=0,nz_plot-1	!z
 iz=z_values_plot(iz_ind)
 WRITE(UNIT=az,FMT='(I2.2)') iz

!all
 OPEN(unit=100,file=trim(label)//'/g/dos_fourier_abs_nz'//az//'_e'//ae2,status='unknown')
! OPEN(unit=101,file=trim(label)//'/g/dos_fourier_img_nz'//az//'_e'//ae2,status='unknown')
! OPEN(unit=102,file=trim(label)//'/g/dos_fourier_rea_nz'//az//'_e'//ae2,status='unknown')
 OPEN(unit=103,file=trim(label)//'/g/dos_fourier_rea_nztot_e'//ae2,status='unknown')

 do ik=0, nkr_plot*nx_kr*ny_kr-1

! if(nx_kr==1 .and. ny_kr==1 .and. mod(nx,2)==1 .and. mod(ny,2)==1) then	!prints real part (qpi is real)
!write(100,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,1),kr_vecs_nscf_qpi(ik,2), &
!				(real(dos_fourier(ik,iz_ind,i)), i=1,10),&	!G-G0
!				(real(dos_fourier(ik,iz_ind,i)), i=21,30)	!G_born-G0
!				real(dos_fourier(ik,iz_ind,1)), real(dos_fourier(ik,iz_ind,2)), real(dos_fourier(ik,iz_ind,3)), real(dos_fourier(ik,iz_ind,4)),&
!				real(dos_fourier(ik,iz_ind,5)), real(dos_fourier(ik,iz_ind,6)),&
!				real(dos_fourier(ik,iz_ind,7)), real(dos_fourier(ik,iz_ind,8)), real(dos_fourier(ik,iz_ind,9))
				
!else	!print modulus (qpi is in general not real)

write(100,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,5),kr_vecs_nscf_qpi(ik,6), (abs(dos_fourier(ik,iz_ind,j)),j=1, entry_dos)
!				(real(dos_fourier(ik,iz_ind,i)), i=1,10),&	!G-G0
!				(real(dos_fourier(ik,iz_ind,i)), i=21,30)	!G_born-G0
				!, abs(dos_fourier(ik,iz_ind,2)), abs(dos_fourier(ik,iz_ind,3)), abs(dos_fourier(ik,iz_ind,4)),&
				!abs(dos_fourier(ik,iz_ind,5)), abs(dos_fourier(ik,iz_ind,6)),&
				!abs(dos_fourier(ik,iz_ind,7)), abs(dos_fourier(ik,iz_ind,8)), abs(dos_fourier(ik,iz_ind,9))
!endif

!if (.not. lfftreal) write(101,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,1),kr_vecs_nscf_qpi(ik,2), (aimag(dos_fourier(ik,iz_ind,j)),j=1, entry_dos)
!				(real(dos_fourier(ik,iz_ind,i)), i=1,10),&	!G-G0
!				(real(dos_fourier(ik,iz_ind,i)), i=21,30)	!G_born-G0
!				aimag(dos_fourier(ik,iz_ind,1)), aimag(dos_fourier(ik,iz_ind,2)), aimag(dos_fourier(ik,iz_ind,3)), aimag(dos_fourier(ik,iz_ind,4)),&
!				aimag(dos_fourier(ik,iz_ind,5)), aimag(dos_fourier(ik,iz_ind,6)),&
!				aimag(dos_fourier(ik,iz_ind,7)), aimag(dos_fourier(ik,iz_ind,8)), aimag(dos_fourier(ik,iz_ind,9))


!if (lfftreal) 
!write(102,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,5),kr_vecs_nscf_qpi(ik,6), (real(dos_fourier(ik,iz_ind,j)),j=1, entry_dos)
!write(103,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,5),kr_vecs_nscf_qpi(ik,6), (sum(real(dos_fourier(ik,:,j))),j=1, entry_dos)
write(103,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,5),kr_vecs_nscf_qpi(ik,6), (real(dos_fourier(ik,0,j))+real(dos_fourier(ik,1,j))/10,j=1, entry_dos)
!				(real(dos_fourier(ik,iz_ind,i)), i=1,10),&	!G-G0
!				(real(dos_fourier(ik,iz_ind,i)), i=21,30)	!G_born-G0
!				real(dos_fourier(ik,iz_ind,1)), real(dos_fourier(ik,iz_ind,2)), real(dos_fourier(ik,iz_ind,3)), real(dos_fourier(ik,iz_ind,4)),&
!				real(dos_fourier(ik,iz_ind,5)), real(dos_fourier(ik,iz_ind,6)),&
!				real(dos_fourier(ik,iz_ind,7)), real(dos_fourier(ik,iz_ind,8)), real(dos_fourier(ik,iz_ind,9))


enddo	!ik
 close(unit=100)
 close(unit=101)
 close(unit=103)


!all_g
if(all_g_fourier) then
 OPEN(unit=100,file=trim(label)//'/g/g_fourier_nz'//az//'_e'//ae2,status='unknown')
! OPEN(unit=101,file=trim(label)//'/g/g_fourier_imag_nz'//az//'_e'//ae2,status='unknown')

 do ik=0, nkr_plot*nx_kr*ny_kr-1
write(100,"(2f10.5,1000e13.5)"), kr_vecs_nscf_qpi(ik,5),kr_vecs_nscf_qpi(ik,6), (real(g_fourier(ik,iz_ind,i)), i=1,2*dim_hk**2)
!write(101,"(2f10.5,1000e13.5)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), (aimag(g_fourier(ik,iz_ind,i)), i=1,2*dim_hk**2)	!should be zero
			
enddo	!ik
 close(unit=100)
! close(unit=101)
endif

enddo	!z

 !summed over z, only the real part (unless it makes no sense)
! OPEN(unit=100,file=trim(label)//'/g/dos_fourier_rea_nztot_e'//ae2,status='unknown')

!do ik=0, nkr_plot*nx_kr*ny_kr-1
!! write(100,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,1),kr_vecs_nscf_qpi(ik,2), (real(sum(dos_fourier(ik,:,j))),j=1,9)	!sum is over z
!if (nz_plot>1) write(100,"(2f10.5,20e14.6)"), kr_vecs_nscf_qpi(ik,1),kr_vecs_nscf_qpi(ik,2), (real(dos_fourier(ik,0,j))+alfa_nz*real(dos_fourier(ik,1,j)),j=1,entry_dos)
!enddo

! close (unit=100)







time_1_w=secnds(0.0)
time_write=time_Write+time_1_w-time_0_w


deallocate(dos_fourier)
if(all_g_fourier) deallocate(g_fourier)

endif !lqpi

deallocate(dos_fin)
if (all_g) deallocate(g_fin)


enddo



if (lqpi)deallocate(kr_vecs_nscf_qpi)
deallocate(V)
!deallocate(w)
!deallocate(e)
!deallocate(w0)
!deallocate(e0)
if (allocated(exp_factor))deallocate(exp_factor)
if (lvkn) deallocate(V_kknn)
deallocate(t_hop)



!write(*,"(a21,f10.2)") "Time V_kn: ",time_vkn
write(*,"(a21,f10.2)") "Time fft: ",time_fft
write(*,"(a21,f10.2)") "Time T: ",time_gg0t
!write(*,"(a21,f10.2)") "Time rhoV**2rho: ",time_rhovrho
write(*,"(a21,f10.2)") "Time G_0(k):",time_g0k
write(*,"(a21,f10.2)") "Time G_0(r-r'):",time_g0r_r
write(*,"(a21,f10.2)") "Time G_0(r,r'):",time_g0rr
write(*,"(a21,f10.2)") "Time G0TG0:",time_g0tg0
write(*,"(a21,f10.2)") "Time DOS:",time_dos
write(*,"(a21,f10.2)") "Time Fourier:",time_fourier
write(*,"(a21,f10.2)") "Time write:",time_write




end subroutine do_green
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Invertcmplx2(N, A, X, B )

INTEGER :: INFO, ITER, N
INTEGER, allocatable :: IPIV( : )
real(kind=idp), allocatable ::  RWORK( : )
COMPLEX (kind=idpc),allocatable::  SWORK( : )
COMPLEX (kind=idpc)::    A( 0:N-1, 0:N-1 ), B( 0:N-1, 0:N-1 ),  X( 0:N-1, 0:N-1 )
COMPLEX (kind=idpc),allocatable	::WORK( :,: ), A2(:,:)

allocate(ipiv(0:N-1))
allocate(rwork(0:N-1))
allocate(swork(0:N*2*N-1))
allocate(work(0:N-1,0:N-1))
allocate(A2(0:N-1,0:N-1))
A2=A

!NRHS=N
!LDA=N
!LDB=N
!LDX=N

time_0_inv=secnds(0.0)

call ZCGESV( N, N, A2, N, IPIV, B, N, X, N, WORK,SWORK, RWORK, ITER, INFO )
!call CGESV( N, N, A, N, IPIV, B, N, X, N, WORK,SWORK, RWORK, ITER, INFO )

time_1_inv=secnds(0.0)

time_invert=time_invert+time_1_inv-time_0_inv


IF (INFO.ne.0) THEN
    WRITE(6,*) ' WARNING: ZCGESV INFO = ', INFO
ENDIF




deallocate(ipiv)
deallocate(rwork)
deallocate(swork)
deallocate(work)
deallocate(A2)


END SUBROUTINE invertcmplx2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Invertreal(N, A, X, B )

INTEGER :: INFO, ITER, N
INTEGER, allocatable :: IPIV( : )
real(kind=idp), allocatable ::  RWORK( : )
real (kind=idp),allocatable::  SWORK( : )
real (kind=idp)::    A( 0:N-1, 0:N-1 ), B( 0:N-1, 0:N-1 ),  X( 0:N-1, 0:N-1 )
real (kind=idp),allocatable	::WORK( :,: )!, A2(:,:)

allocate(ipiv(0:N-1))
allocate(rwork(0:N-1))
allocate(swork(0:N*2*N-1))
allocate(work(0:N-1,0:N-1))
!allocate(A2(0:N-1,0:N-1))
!A2=A

!NRHS=N
!LDA=N
!LDB=N
!LDX=N

time_0_inv=secnds(0.0)

!call ZCGESV( N, N, A2, N, IPIV, B, N, X, N, WORK,SWORK, RWORK, ITER, INFO )
call CGESV( N, N, A, N, IPIV, B, N, X, N, WORK,SWORK, RWORK, ITER, INFO )

time_1_inv=secnds(0.0)

time_invert=time_invert+time_1_inv-time_0_inv


IF (INFO.ne.0) THEN
    WRITE(6,*) ' WARNING: ZCGESV INFO = ', INFO
ENDIF




deallocate(ipiv)
deallocate(rwork)
deallocate(swork)
deallocate(work)
!deallocate(A2)


END SUBROUTINE invertreal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE Invertcmplx(N, A, X, B )

INTEGER :: INFO, ITER, N
INTEGER, allocatable :: IPIV( : )
real(kind=idp), allocatable ::  RWORK( : )
COMPLEX (kind=idpc),allocatable::  SWORK( :,: )
COMPLEX (kind=idpc)::    A( 0:N-1, 0:N-1 ), B( 0:N-1, 0:N-1 ),  X( 0:N-1, 0:N-1 )
COMPLEX (kind=idpc),allocatable	::WORK( :,: ), A2(:,:), B2(:,:), X2(:,:)

allocate(ipiv(N))
allocate(rwork(N))
allocate(swork(N,2*N))
allocate(work(N,N))
allocate(A2(N,N))
allocate(B2(N,N))
allocate(X2(N,N))

do i=0,N-1
do j=0,N-1
a2(i+1,j+1)=a(i,j)
b2(i+1,j+1)=b(i,j)
enddo
enddo


!NRHS=N
!LDA=N
!LDB=N
!LDX=N

time_0_inv=secnds(0.0)

call ZCGESV( N, N, A2, N, IPIV, B2, N, X2, N, WORK,SWORK, RWORK, ITER, INFO )
!call CGESV( N, N, A, N, IPIV, B, N, X, N, WORK,SWORK, RWORK, ITER, INFO )

time_1_inv=secnds(0.0)

time_invert=time_invert+time_1_inv-time_0_inv


IF (INFO.ne.0) THEN
    WRITE(6,*) ' WARNING: ZCGEV INFO = ', INFO
ENDIF


do i=0,N-1
do j=0,N-1
x(i,j)=x2(i+1,j+1)
enddo
enddo




deallocate(ipiv)
deallocate(rwork)
deallocate(swork)
deallocate(work)
deallocate(A2)
deallocate(B2)
deallocate(X2)


END SUBROUTINE invertcmplx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matmul1(N, A, B, X)
COMPLEX (kind=idpc)::    A( 0:N-1, 0:N-1 ), B( 0:N-1, 0:N-1 ),  X( 0:N-1, 0:N-1 )
COMPLEX (kind=idpc),allocatable::    A2(:,:), B2(:,:),  X2(:,:)
integer	::N,i,j


!allocate(A2(N,N))
!allocate(B2(N,N))
!allocate(X2(N,N))
!A2=0
!B2=0
!X2=0


do i=0,N-1
do j=0,N-1
!A2(i+1,j+1)=A(i,j)
!B2(i+1,j+1)=B(i,j)
enddo
enddo

time_0_matmul=secnds(0.0)

!CALL ZGEMM( 'N', 'N', N,N, N, 1.0D0, A2, N, B2, N, 0.0D0, X2, N )
CALL ZGEMM( 'N', 'N', N,N, N, 1.0D0, A, N, B, N, 0.0D0, X, N )
time_1_matmul=secnds(0.0)

time_matmul=time_matmul+time_1_matmul-time_0_matmul


do i=0,N-1
do j=0,N-1
!X(i,j)=X2(i+1,j+1)
enddo
enddo



!deallocate(A2)
!deallocate(B2)
!deallocate(X2)


end subroutine matmul1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matmul_hc(N, A, B, X)
COMPLEX (kind=idpc)::    A( 0:N-1, 0:N-1 ), B( 0:N-1, 0:N-1 ),  X( 0:N-1, 0:N-1 )
COMPLEX (kind=idpc),allocatable::    A2(:,:), B2(:,:),  X2(:,:)
integer	::N

time_0_matmul=secnds(0.0)

CALL ZGEMM( 'N', 'C', N,N, N, 1.0D0, A, N, B, N, 0.0D0, X, N )
time_1_matmul=secnds(0.0)

time_matmul=time_matmul+time_1_matmul-time_0_matmul

end subroutine matmul_hc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!A(M,K)*B(K,N)=X(M,N)
subroutine matmul2(M,N,K, A, B, X)
COMPLEX (kind=idpc)::    A( 0:M-1, 0:K-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
integer	::M,N,K,i,j
COMPLEX (kind=idpc),allocatable::    A2(:,:), B2(:,:),  X2(:,:)



!allocate(A2(M,K))
!allocate(B2(K,N))
!allocate(X2(M,N))
!A2=0
!B2=0
!X2=0


!do i=0,M-1
!do j=0,K-1
!A2(i+1,j+1)=A(i,j)
!enddo
!enddo
!do i=0,K-1
!do j=0,N-1
!B2(i+1,j+1)=B(i,j)
!enddo
!enddo



time_0_matmul=secnds(0.0)

!CALL ZGEMM( 'N', 'N', M,N, K, 1.0D0, A2, M, B2, K, 0.0D0, X2, M )
CALL ZGEMM( 'N', 'N', M,N, K, 1.0D0, A, M, B, K, 0.0D0, X, M )
time_1_matmul=secnds(0.0)

time_matmul=time_matmul+time_1_matmul-time_0_matmul


!do i=0,M-1
!do j=0,N-1
!X(i,j)=X2(i+1,j+1)
!enddo
!enddo


!deallocate(A2)
!deallocate(B2)
!deallocate(X2)



end subroutine matmul2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!A(M,K)*B(K,N)=X(M,N)
subroutine matmul2hc(M,N,K, A, B, X)
COMPLEX (kind=idpc)::    A( 0:M-1, 0:K-1 ), B( 0:N-1, 0:K-1 ),  X( 0:M-1, 0:N-1 )
integer	::M,N,K,i,j
COMPLEX (kind=idpc),allocatable::    A2(:,:), B2(:,:),  X2(:,:)



!allocate(A2(M,K))
!allocate(B2(K,N))
!allocate(X2(M,N))
!A2=0
!B2=0
!X2=0


!do i=0,M-1
!do j=0,K-1
!A2(i+1,j+1)=A(i,j)
!enddo
!enddo
!do i=0,K-1
!do j=0,N-1
!B2(i+1,j+1)=B(i,j)
!enddo
!enddo



time_0_matmul=secnds(0.0)

!CALL ZGEMM( 'N', 'N', M,N, K, 1.0D0, A2, M, B2, K, 0.0D0, X2, M )
CALL ZGEMM( 'N', 'C', M,N, K, 1.0D0, A, M, B, N, 0.0D0, X, M )
time_1_matmul=secnds(0.0)

time_matmul=time_matmul+time_1_matmul-time_0_matmul


!do i=0,M-1
!do j=0,N-1
!X(i,j)=X2(i+1,j+1)
!enddo
!enddo


!deallocate(A2)
!deallocate(B2)
!deallocate(X2)



end subroutine matmul2hc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!conjgTr(A(K,M))*B(K,N)=X(M,N) 
subroutine matmul2hc2(M,N,K, A, B, X)
COMPLEX (kind=idpc)::    A( 0:K-1, 0:M-1 ), B( 0:K-1, 0:N-1 ),  X( 0:M-1, 0:N-1 )
integer	::M,N,K,i,j
COMPLEX (kind=idpc),allocatable::    A2(:,:), B2(:,:),  X2(:,:)



!allocate(A2(M,K))
!allocate(B2(K,N))
!allocate(X2(M,N))
!A2=0
!B2=0
!X2=0


!do i=0,M-1
!do j=0,K-1
!A2(i+1,j+1)=A(i,j)
!enddo
!enddo
!do i=0,K-1
!do j=0,N-1
!B2(i+1,j+1)=B(i,j)
!enddo
!enddo



time_0_matmul=secnds(0.0)

!CALL ZGEMM( 'N', 'N', M,N, K, 1.0D0, A2, M, B2, K, 0.0D0, X2, M )
CALL ZGEMM( 'C', 'N', M,N, K, 1.0D0, A, K, B, K, 0.0D0, X, M )
time_1_matmul=secnds(0.0)

time_matmul=time_matmul+time_1_matmul-time_0_matmul


!do i=0,M-1
!do j=0,N-1
!X(i,j)=X2(i+1,j+1)
!enddo
!enddo


!deallocate(A2)
!deallocate(B2)
!deallocate(X2)



end subroutine matmul2hc2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kvecs_r_arpes(k_vecs_r_arpes,nk_r_arpes,list_k_vecs_r_arpes)	!nscf for arpes
real(kind=idp), allocatable,intent(out) :: k_vecs_r_arpes(:,:)
integer, intent(out)	:: nk_r_arpes
integer, allocatable,intent(out) :: list_k_vecs_r_arpes(:,:)
real(kind=idp)	::k_vec_temp(2)
logical :: ok

!works bor both nx and nkx_r_nscf odd


print *, "Setting k points for arpes"

nk_r_arpes=(nx*nkx_r_nscf)**2

allocate(k_vecs_r_arpes(0:nk_r_arpes-1,4))
allocate(list_k_vecs_r_arpes(0:nk_r_arpes-1,6))

k_vecs_r_arpes=-100
list_k_vecs_r_arpes=-100

!goto 100

nkx_l=nx/2
nky_l=ny/2

if(lysameasx) ysameasx=1
if( .not. lysameasx) ysameasx=0

ik_arpes=-1

do ik=0, nk_r-1
do sy=1,-1,-2		!y->-y
 do sx=1,-1,-2			!x->-x
  do swap_xy=0,ysameasx		!x->y
   do ikx=-nkx_l,nkx_l			!Kx large
    do iky=-nkx_l,nky_l		!Ky large
    
     if (swap_xy==0) k_vec_temp(1)=sx*(k_vecs_r(ik,1)+ikx)/nx		
     if (swap_xy==0) k_vec_temp(2)=sy*(k_vecs_r(ik,2)+iky)/nx

     if (swap_xy==1) k_vec_temp(1)=sx*(k_vecs_r(ik,2)+iky)/nx		
     if (swap_xy==1) k_vec_temp(2)=sy*(k_vecs_r(ik,1)+ikx)/nx

!checks if k point already exists
 ok=.true.
 do j=0, ik_arpes
 if(abs(k_vec_temp(1)-k_vecs_r_arpes(j,1))+abs(k_vec_temp(2)-k_vecs_r_arpes(j,2))<1e-5) ok=.false.
 enddo

if (ok) then
ik_arpes=ik_arpes+1
print *, ik_arpes
k_vecs_r_arpes(ik_arpes,1)=k_vec_temp(1)
k_vecs_r_arpes(ik_arpes,2)=k_vec_temp(2)
list_k_vecs_r_arpes(ik_arpes,1)=ik
list_k_vecs_r_arpes(ik_arpes,2)=ikx
list_k_vecs_r_arpes(ik_arpes,3)=iky
list_k_vecs_r_arpes(ik_arpes,4)=sx
list_k_vecs_r_arpes(ik_arpes,5)=sy
list_k_vecs_r_arpes(ik_arpes,6)=swap_xy
endif

enddo
enddo
enddo
enddo
enddo
enddo

print *, nk_r_arpes, ik_arpes

write(*, "(7a4,6a10)"), "i", "ik","ikx", "iky", "sx", "sy", "s_xy", "k_x", "k_y", "K_x", "K_y", "Ktot_x", "Ktot_y"
do i=0, nk_r_arpes-1
write(*, "(7i4,6f10.5)"), i, list_k_vecs_r_arpes(i,1),list_k_vecs_r_arpes(i,2),list_k_vecs_r_arpes(i,3),list_k_vecs_r_arpes(i,4),list_k_vecs_r_arpes(i,5),list_k_vecs_r_arpes(i,6), &
k_vecs_r(list_k_vecs_r_arpes(i,1),1),k_vecs_r(list_k_vecs_r_arpes(i,1),2),real(list_k_vecs_r_arpes(i,2))/nx,real(list_k_vecs_r_arpes(i,2))/ny,  k_vecs_r_arpes(i,1),k_vecs_r_arpes(i,2)
enddo

!stop
!100 continue


end subroutine set_kvecs_r_arpes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_g0xtmp(x,y,n,m)
integer	::x,y,n,m
integer	::indx,indy

!x=ix-ix2
!y=iy-iy2
!n=nkx_plot
!m=nky_plot

index_g0xtmp=mod(x+n,n)+mod(y+m,m)*n

!!!!!!!!!!!!!!!!

!if (mod(abs(x),n)<m) indx=mod(abs(x), n)
!if (mod(abs(x),n)>m-1) indx=n-mod(abs(x),n)
!if (mod(abs(y),n)<m) indy=mod(abs(y), n)
!if (mod(abs(y),n)>m-1) indy=n-mod(abs(y),n)

!index_g0xtmp=indx+m*indy

!print *, n,m,x,y,indx,indy

end function index_g0xtmp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fortran example.
! 2D complex to complex, and real to conjugate-even
subroutine fft
!Use MKL_DFTI
Complex ::  X_2D(32,100)
Real :: Y_2D(34, 102)
Complex ::  X(3200)
Real :: Y(3468)
Equivalence (X_2D, X)
Equivalence (Y_2D, Y)
!type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
Integer :: Status, L(2)
!...put input data into X_2D(j,k), Y_2D(j,k), 1<=j=32,1<=k<=100
!...set L(1) = 32, L(2) = 100
!...the transform is a 32-by-100
 
! Perform a complex to complex transform
!Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE,&
!          DFTI_COMPLEX, 2, L)
!Status = DftiCommitDescriptor( My_Desc1_Handle)
!Status = DftiComputeForward( My_Desc1_Handle, X)
!Status = DftiFreeDescriptor(My_Desc1_Handle)
! result is given by X_2D(j,k), 1<=j<=32, 1<=k<=100
 
! Perform a real to complex conjugate-even transform
!Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_SINGLE,&
!          DFTI_REAL, 2, L)
!Status = DftiCommitDescriptor( My_Desc2_Handle)
!Status = DftiComputeForward( My_Desc2_Handle, Y)
!Status = DftiFreeDescriptor(My_Desc2_Handle)
! result is given by the complex value z(j,k) 1<=j<=32; 1<=k<=100
! and is stored in CCS format!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine index_kkq(ik,ik2,iq)
integer, intent(in)	::ik,ik2
integer, intent(out)	::iq
real::kx,ky
integer::i,ikx,iky,ikx2,iky2,iqx,iqy

ikx=mod(ik,nkx_plot)
iky=(ik-ikx)/nkx_plot
ikx2=mod(ik2,nkx_plot)
iky2=(ik2-ikx2)/nkx_plot

iqx=mod(ikx-ikx2+3*nkx_plot/2,nkx_plot)
iqy=mod(iky-iky2+3*nkx_plot/2,nkx_plot)

iq=iqx+iqy*nkx_plot

!if (abs(kvecs(ik,1)- kvecs(ik2,1)-kvecs(iq,1))>1e-5 .and. abs(kvecs(ik,1)- kvecs(ik2,1)-kvecs(iq,1)-1)>1e-4 .and. abs(kvecs(ik,1)- kvecs(ik2,1)-kvecs(iq,1))>1e-5) then
!write(*, "(a3,4f10.5,5i5)") "x",kr_vecs_nscf(ik,1),kr_vecs_nscf(ik2,1),kr_vecs_nscf(ik,1)- kr_vecs_nscf(ik2,1),kr_vecs_nscf(iq,1), ikx,ikx2,iqx
!endif

!if (abs(kvecs(ik,2)- kvecs(ik2,2)-kvecs(iq,2))>1e-5 .and. abs(kvecs(ik,2)- kvecs(ik2,2)-kvecs(iq,2)-1)>1e-4 .and. abs(kvecs(ik,2)- kvecs(ik2,2)-kvecs(iq,2))>1e-5) then
!write(*, "(a3,4f10.5,5i5)") "y",kr_vecs_nscf(ik,2),kr_vecs_nscf(ik2,2),kr_vecs_nscf(ik,2)- kr_vecs_nscf(ik2,2),kr_vecs_nscf(iq,2), iky,iky2,iqy
!endif


end subroutine index_kkq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine w_kr_analysis(e_kr_tot_nscf,w_kr_tot_nscf)
complex(KIND=idpc),intent(in) 	:: w_kr_tot_nscf(0:size_h_kr_red-1,0:size_h_kr_red_e-1,0:nkr_plot-1)
real (kind=idp),intent(in)	:: e_kr_tot_nscf(0:size_h_kr-1,0:nkr_plot-1)
integer	::ikx,iky,ik,i,i2,ikx2,iky2,ik2
real (kind=idp), allocatable	:: sigma_xyz(:,:,:,:,:)
!				 is  is2  icf xyz
integer	::is,is2,icf,ixyz
 CHARACTER(2) :: az
 CHARACTER(5) :: ak

 print *, "w_kr_analysis"


!			k		l   	z    c,f,cf  xyz
allocate(sigma_xyz(0:nkr_plot-1, 0:size_h_kr-1, 0:nz_scf-1, 0:1, 1:3)) 
sigma_xyz=0

do ik=0, nkr_plot-1
 do i=0,size_h_kr_red_e-1
  do iz=0, nz_scf-1
   do icf=0,1
    do is=0,1
     do is2=0,1
      do ixyz=1,3
       do o1=0,l_orbital-1
        do o2=0,l_orbital-1

ind= iz + icf*nz_scf + is*2*nz_scf  + o1*4*nz_scf
ind2=iz + icf*nz_scf + is2*2*nz_scf + o2*4*nz_scf

sigma_xyz(ik,i,iz,icf,ixyz)=	sigma_xyz(ik,i,iz,icf,ixyz)+&
				conjg(w_kr_tot_nscf(ind,i,ik))*w_kr_tot_nscf(ind2,i,ik)*sigma_xyz_cf(is,is2,o1,o2,icf,ixyz)
!				abs(w_kr_tot_nscf(ind,i,ik))**2*kr_vecs_nscf(ik,4)/nkx_plot/nky_plot
   
   
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo



!z resolved (central bands)
do iz=0, nz_scf-1
 WRITE(UNIT=az,FMT='(I2.2)') iz

OPEN(unit=300,file=trim(label)//'/sf/spin_m_z'//az,status='unknown')
OPEN(unit=301,file=trim(label)//'/sf/spin_p_z'//az,status='unknown')

do ik=0, nkr_plot-1

write(300,"(2f10.5,12e16.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), &
				sigma_xyz(ik,2*nz-2,iz,0,1)+sigma_xyz(ik,2*nz-1,iz,0,1),sigma_xyz(ik,2*nz-2,iz,0,2)+sigma_xyz(ik,2*nz-1,iz,0,2),sigma_xyz(ik,2*nz-2,iz,0,3)+sigma_xyz(ik,2*nz-1,iz,0,3),&
				sigma_xyz(ik,2*nz-2,iz,1,1)+sigma_xyz(ik,2*nz-1,iz,1,1),sigma_xyz(ik,2*nz-2,iz,1,2)+sigma_xyz(ik,2*nz-1,iz,1,2),sigma_xyz(ik,2*nz-2,iz,1,3)+sigma_xyz(ik,2*nz-1,iz,1,3)
				
write(301,"(2f10.5,12e16.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), &
				sigma_xyz(ik,2*nz,iz,0,1)+sigma_xyz(ik,2*nz+1,iz,0,1),sigma_xyz(ik,2*nz,iz,0,2)+sigma_xyz(ik,2*nz+1,iz,0,2),sigma_xyz(ik,2*nz,iz,0,3)+sigma_xyz(ik,2*nz+1,iz,0,3),&
				sigma_xyz(ik,2*nz,iz,1,1)+sigma_xyz(ik,2*nz+1,iz,1,1),sigma_xyz(ik,2*nz,iz,1,2)+sigma_xyz(ik,2*nz+1,iz,1,2),sigma_xyz(ik,2*nz,iz,1,3)+sigma_xyz(ik,2*nz+1,iz,1,3)

enddo
 close (unit=300)
 close (unit=301)
					
enddo


!summed over z

OPEN(unit=300,file=trim(label)//'/sf/spin_m',status='unknown')
OPEN(unit=301,file=trim(label)//'/sf/spin_p',status='unknown')


do ik=0, nkr_plot-1
write(300,"(2f10.5,12e16.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), &
	sum(sigma_xyz(ik,2*nz-2,:,0,1))+sum(sigma_xyz(ik,2*nz-1,:,0,1)),sum(sigma_xyz(ik,2*nz-2,:,0,2))+sum(sigma_xyz(ik,2*nz-1,:,0,2)),sum(sigma_xyz(ik,2*nz-2,:,0,3))+sum(sigma_xyz(ik,2*nz-1,:,0,3)),&
	sum(sigma_xyz(ik,2*nz-2,:,1,1))+sum(sigma_xyz(ik,2*nz-1,:,1,1)),sum(sigma_xyz(ik,2*nz-2,:,1,2))+sum(sigma_xyz(ik,2*nz-1,:,1,2)),sum(sigma_xyz(ik,2*nz-2,:,1,3))+sum(sigma_xyz(ik,2*nz-1,:,1,3))
				
write(301,"(2f10.5,12e16.8)"), kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2), &
	sum(sigma_xyz(ik,2*nz,:,0,1))+sum(sigma_xyz(ik,2*nz+1,:,0,1)),sum(sigma_xyz(ik,2*nz,:,0,2))+sum(sigma_xyz(ik,2*nz+1,:,0,2)),sum(sigma_xyz(ik,2*nz,:,0,3))+sum(sigma_xyz(ik,2*nz+1,:,0,3)),&
	sum(sigma_xyz(ik,2*nz,:,1,1))+sum(sigma_xyz(ik,2*nz+1,:,1,1)),sum(sigma_xyz(ik,2*nz,:,1,2))+sum(sigma_xyz(ik,2*nz+1,:,1,2)),sum(sigma_xyz(ik,2*nz,:,1,3))+sum(sigma_xyz(ik,2*nz+1,:,1,3))
enddo

 close (unit=300)
 close (unit=301)

!stop



!DOS

  do iz_ind=0,nz_plot-1	!z
   iz=z_values_plot(iz_ind)
 do ik=0, nkr_plot-1			!kpoint

 WRITE(UNIT=ak,FMT='(I5.5)') ik
 WRITE(UNIT=az,FMT='(I2.2)') iz
 
 OPEN(unit=100,file=trim(label)//'/sf/sf_kr_spins_nz'//az//'_k'//ak,status='unknown')

  
  write(100,"(i10)") size_h_kr

  do i=0, size_h_kr-1		!eigenvalue 
    write(100,format_smear) e_kr_tot_nscf(i,ik), &
    				sigma_xyz(ik,i,iz,0,1),sigma_xyz(ik,i,iz,0,2),sigma_xyz(ik,i,iz,0,3),&	!cx, cy, cz
				sigma_xyz(ik,i,iz,1,1),sigma_xyz(ik,i,iz,1,2),sigma_xyz(ik,i,iz,1,3)	!fx, fy, fz

  enddo
    write(100,format_smear) kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2),kr_vecs_nscf(ik,4)

 close(unit=100)


 enddo
enddo












deallocate(sigma_xyz)

!ikx=nkx_plot/2-3
!iky=0
!ik=ikx+nkx_plot*iky
!ikx2=nkx_plot/2+3
!iky2=nkx_plot-1
!ik2=ikx2+nkx_plot*iky2

!write(*, "(2i6,2f10.5, 4f12.6)"),ikx,iky,kr_vecs_nscf(ik,1),kr_vecs_nscf(ik,2)
!write(*, "(2i6,2f10.5, 4f12.6)"),ikx2,iky2,kr_vecs_nscf(ik2,1),kr_vecs_nscf(ik2,2)


!do i=0,1		
! i2=i+2*nz
! write (*, "(2i5,10f15.10)"), i,i2,e_kr_tot_nscf(i2,ik),e_kr_tot_nscf(i2,ik2)

!   do iz=0, nz-1
!    do icf=0,1
!     do is=0,1
!ind=iz + icf*nz + is*2*nz
!if (iz<2) write (*, "(3i5,10f15.10)"),iz,icf,is, real(w_kr_tot_nscf(ind,i2,ik)), aimag(w_kr_tot_nscf(ind,i2,ik)),real(w_kr_tot_nscf(ind,i2,ik2)), aimag(w_kr_tot_nscf(ind,i2,ik2))
!enddo
!enddo
!enddo
!enddo



end subroutine w_kr_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function phase_of_exp(x)
!
! given x = exp(2*i*phi), finds phi
!
!  USE kinds, only : DP

  implicit none

  real(KIND=idp) :: phase_of_exp, delta
  complex (KIND=idpc):: x

  delta = DBLE(x)
  if(delta.gt.1.d0) delta = 1.d0
  if(delta.lt.-1.d0) delta = -1.d0
  delta = DACOS(delta)

  IF(AIMAG(x).lt.0.d0) delta = -delta
  phase_of_exp = 0.5d0 * delta

  return
end function phase_of_exp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function phase_of_t(x)

  implicit none

  real(KIND=idp) :: phase_of_t
  complex (KIND=idpc):: x

!  delta = DBLE(x)
! if (real(x)>=0) then 
! phase_of_t=atan(aimag(x)/real(x))
 phase_of_t=atan2(aimag(x),real(x))
! else if (real(x)<0 .and. aimag(x)>=0) then
! phase_of_t=atan(aimag(x)/real(x))+pi
! else
! phase_of_t=atan(aimag(x)/real(x))-pi
! endif
 

  return
end function phase_of_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_spin
real (kind=idp):: c2t, s2t
integer	:: ixyz

allocate(sigma_xyz_cf(0:1,0:1,0:l_orbital_f-1,0:l_orbital_f-1,0:1, 1:3))
allocate(sigma_pauli(0:1,0:1,1:3))


!0:down
!1:up

sigma_pauli=0
!sigma_x
sigma_pauli(1,0,1)=1
sigma_pauli(0,1,1)=1
!sigma_y
sigma_pauli(1,0,2)=-(0,1)
sigma_pauli(0,1,2)=+(0,1)
!sigma_z
sigma_pauli(0,0,3)=-1
sigma_pauli(1,1,3)=1


sigma_xyz_cf=0


!f
do o1=0, l_orbital_f-1
do o2=0, l_orbital_f-1
do is=0,1	!pseudospin	
do is2=0,1	!pseudospin
do ixyz=1,3
do icf=1,1
do k=-3,3	!m, to be traced
do is3=0,1	!spin
do is4=0,1	!spin

sigma_xyz_cf(is,is2,o1,o2,icf,ixyz)=sigma_xyz_cf(is,is2,o1,o2,icf,ixyz)+alpha_msigma(o1,is,is3,k)*sigma_pauli(is3,is4,ixyz)*alpha_msigma(o2,is2,is4,k)

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo

!s-d
do o1=0, l_orbital-1
do is=0,1	!pseudospin=spin	
do is2=0,1	!pseudospin=spin
do ixyz=1,3
do icf=0,0

sigma_xyz_cf(is,is2,o1,o1,icf,ixyz)=sigma_pauli(is,is2,ixyz)
!sigma_xyz_cf(is,is2,o1,o1,1,ixyz)=sigma_pauli(is,is2,ixyz)

enddo
enddo
enddo
enddo
enddo


!print 
!s-d


if (l2dmodel .or. phi_type=="q2d") then	!here I use the simples spin

sigma_xyz_cf=0

do is=0,1	!pseudospin=spin	
do is2=0,1	!pseudospin=spin
do ixyz=1,3

sigma_xyz_cf(is,is2,0,0,1,ixyz)=sigma_pauli(is,is2,ixyz)

enddo
enddo
enddo
endif








print *, "Spin on Basis:"
print *, "s-d"
icf=0

if (l_orbital==2) then

ixyz=1
write (*, "(a5, 4a18)"), "s_x", "d1-", "d1+", "d2-", "d2+"
write (*, "(a5, 4(f10.4,f8.4))"), "d1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)

ixyz=2
write (*, "(a5, 4a18)"), "s_y", "d1-", "d1+", "d2-", "d2+"
write (*, "(a5, 4(f10.4,f8.4))"), "d1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)

ixyz=3
write (*, "(a5, 4a18)"), "s_z", "d1-", "d1+", "d2-", "d2+"
write (*, "(a5, 4(f10.4,f8.4))"), "d1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)

print *, ""

else

ixyz=1
write (*, "(a5, 4a18)"), "s_x", "d-", "d+"
write (*, "(a5, 4(f10.4,f8.4))"), "d-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)

ixyz=2
write (*, "(a5, 4a18)"), "s_y", "d-", "d+"
write (*, "(a5, 4(f10.4,f8.4))"), "d-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)

ixyz=3
write (*, "(a5, 4a18)"), "s_z", "d-", "d+"
write (*, "(a5, 4(f10.4,f8.4))"), "d-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "d+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)

print *, ""

endif


print *, "f"
icf=1

if (l_orbital_f==2) then

ixyz=1
write (*, "(a5, 4a18)"), "s_x", "f1-", "f1+", "f2-", "f2+"
write (*, "(a5, 4(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)

ixyz=2
write (*, "(a5, 4a18)"), "s_y", "d1-", "d1+", "d2-", "d2+"
write (*, "(a5, 4(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)

ixyz=3
write (*, "(a5, 4a18)"), "s_z", "f1-", "f1+", "f2-", "f2+"
write (*, "(a5, 4(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz)


else if (l_orbital_f==3) then
ixyz=1
write (*, "(a5, 6a18)"), "s_x", "f1-", "f1+", "f2-", "f2+", "f3-", "f3+"
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,0,2,icf,ixyz),sigma_xyz_cf(0,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,0,2,icf,ixyz),sigma_xyz_cf(1,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,1,2,icf,ixyz),sigma_xyz_cf(0,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,1,2,icf,ixyz),sigma_xyz_cf(1,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,2,0,icf,ixyz),sigma_xyz_cf(0,1,2,0,icf,ixyz),sigma_xyz_cf(0,0,2,1,icf,ixyz),sigma_xyz_cf(0,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,2,2,icf,ixyz),sigma_xyz_cf(0,1,2,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,2,0,icf,ixyz),sigma_xyz_cf(1,1,2,0,icf,ixyz),sigma_xyz_cf(1,0,2,1,icf,ixyz),sigma_xyz_cf(1,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,2,2,icf,ixyz),sigma_xyz_cf(1,1,2,2,icf,ixyz)

ixyz=2
write (*, "(a5, 4a18)"), "s_y", "f1-", "f1+", "f2-", "f2+","f3-", "f3+"
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,0,2,icf,ixyz),sigma_xyz_cf(0,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,0,2,icf,ixyz),sigma_xyz_cf(1,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,1,2,icf,ixyz),sigma_xyz_cf(0,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,1,2,icf,ixyz),sigma_xyz_cf(1,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,2,0,icf,ixyz),sigma_xyz_cf(0,1,2,0,icf,ixyz),sigma_xyz_cf(0,0,2,1,icf,ixyz),sigma_xyz_cf(0,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,2,2,icf,ixyz),sigma_xyz_cf(0,1,2,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,2,0,icf,ixyz),sigma_xyz_cf(1,1,2,0,icf,ixyz),sigma_xyz_cf(1,0,2,1,icf,ixyz),sigma_xyz_cf(1,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,2,2,icf,ixyz),sigma_xyz_cf(1,1,2,2,icf,ixyz)
ixyz=3
write (*, "(a5, 4a18)"), "s_z", "f1-", "f1+", "f2-", "f2+", "f3-", "f3+"
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz),sigma_xyz_cf(0,0,0,1,icf,ixyz),sigma_xyz_cf(0,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,0,2,icf,ixyz),sigma_xyz_cf(0,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz),sigma_xyz_cf(1,0,0,1,icf,ixyz),sigma_xyz_cf(1,1,0,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,0,2,icf,ixyz),sigma_xyz_cf(1,1,0,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2-", sigma_xyz_cf(0,0,1,0,icf,ixyz),sigma_xyz_cf(0,1,1,0,icf,ixyz),sigma_xyz_cf(0,0,1,1,icf,ixyz),sigma_xyz_cf(0,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,1,2,icf,ixyz),sigma_xyz_cf(0,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f2+", sigma_xyz_cf(1,0,1,0,icf,ixyz),sigma_xyz_cf(1,1,1,0,icf,ixyz),sigma_xyz_cf(1,0,1,1,icf,ixyz),sigma_xyz_cf(1,1,1,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,1,2,icf,ixyz),sigma_xyz_cf(1,1,1,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1-", sigma_xyz_cf(0,0,2,0,icf,ixyz),sigma_xyz_cf(0,1,2,0,icf,ixyz),sigma_xyz_cf(0,0,2,1,icf,ixyz),sigma_xyz_cf(0,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(0,0,2,2,icf,ixyz),sigma_xyz_cf(0,1,2,2,icf,ixyz)
write (*, "(a5, 6(f10.4,f8.4))"), "f1+", sigma_xyz_cf(1,0,2,0,icf,ixyz),sigma_xyz_cf(1,1,2,0,icf,ixyz),sigma_xyz_cf(1,0,2,1,icf,ixyz),sigma_xyz_cf(1,1,2,1,icf,ixyz),&
					 sigma_xyz_cf(1,0,2,2,icf,ixyz),sigma_xyz_cf(1,1,2,2,icf,ixyz)



else

ixyz=1
write (*, "(a5, 4a18)"), "s_x", "f-", "f+"
write (*, "(a5, 4(f10.4,f8.4))"), "f-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)


ixyz=2
write (*, "(a5, 4a18)"), "s_y", "f-", "f+"
write (*, "(a5, 4(f10.4,f8.4))"), "f-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)

ixyz=3
write (*, "(a5, 4a18)"), "s_z", "f-", "f+"
write (*, "(a5, 4(f10.4,f8.4))"), "f-", sigma_xyz_cf(0,0,0,0,icf,ixyz),sigma_xyz_cf(0,1,0,0,icf,ixyz)
write (*, "(a5, 4(f10.4,f8.4))"), "f+", sigma_xyz_cf(1,0,0,0,icf,ixyz),sigma_xyz_cf(1,1,0,0,icf,ixyz)

endif

print *, "*******"


!c,x
!sigma_xyz_cf(0,1,0,1)=1
!sigma_xyz_cf(1,0,0,1)=1
!c,y
!sigma_xyz_cf(0,1,0,2)=+(0,1)	!<down|sigma|up>
!sigma_xyz_cf(1,0,0,2)=-(0,1)	!<up|sigma|down>
!c,z
!sigma_xyz_cf(0,0,0,3)=-1	!down
!sigma_xyz_cf(1,1,0,3)=1		!up


!if (phi_type=="1-2") then

! c2t=3./7.
! s2t=4./7.

!same as above, but with -c2t for x,y, and c2t-s2t for z
!f,x
!sigma_xyz_cf(0,1,1,1)=-c2t
!sigma_xyz_cf(1,0,1,1)=-c2t
!f,y
!sigma_xyz_cf(0,1,1,2)=-(0,1)*c2t	!<-|sigma|+>
!sigma_xyz_cf(1,0,1,2)=+(0,1)*c2t	!<+|sigma|->
!f,z
!sigma_xyz_cf(0,0,1,3)=s2t-c2t
!sigma_xyz_cf(1,1,1,3)=c2t-s2t


!elseif (phi_type=="srb") then

! c2t=2./3.
! s2t=1./3.

!p,x
!sigma_xyz_cf(0,1,1,1)=-c2t
!sigma_xyz_cf(1,0,1,1)=-c2t
!p,y
!sigma_xyz_cf(0,1,1,2)=-(0,1)*c2t	!<-|sigma|+>
!sigma_xyz_cf(1,0,1,2)=+(0,1)*c2t	!<+|sigma|->
!p,z
!sigma_xyz_cf(0,0,1,3)=s2t-c2t
!sigma_xyz_cf(1,1,1,3)=c2t-s2t


!endif


!print *, sigma_xyz_cf


!reflection symmetry operator
n_sym_arpes=13
allocate(simm_refl(0:1,0:1,0:l_orbital_f-1,0:l_orbital_f-1,0:1,0:n_sym_arpes))
allocate(simm_spin(0:1,0:1,0:n_sym_arpes))
allocate(simm_orb(-3:3,-3:3,0:n_sym_arpes))
allocate(simm_orb_d(0:l_orbital-1,0:l_orbital-1,0:n_sym_arpes))
allocate(simm_oper(0:dim_hk-1,0:dim_hk-1,0:n_sym_arpes))
allocate(simm_oper2(0:dim_hk-1,0:dim_hk-1,0:n_sym_arpes))

simm_refl=0
simm_spin=0
simm_orb=0
simm_orb_d=0


print *, "0: Identity"
print *, "1: Parity"
print *, "2: TRS"
print *, "3: M_yz/i"
print *, "4: M_xz/i"
print *, "5: M_xy/i"
print *, "6: M_xz-yz/i"
print *, "7: M_xz+yz/i"
print *, "8: R_z pi/2"
print *, "9: R_z pi"
print *,"10: R_x pi"
print *,"11: R_(x-y) pi"
print *,"12: M_(yx-zx)"

!print *, "8: M_xy-zy"
!print *, "9: M_xy+zy"
!print *, "10: M_yx-zx"
!print *, "11: M_yx+zx"

!0: Identity
j=0
simm_spin(0,0,j)=1
simm_spin(1,1,j)=1
do i=-3,3
simm_orb(i,i,j)=1
enddo
do i=0,l_orbital-1
simm_orb_d(i,i,j)=1
enddo
!1: Parity
j=1
simm_spin(0,0,j)=1
simm_spin(1,1,j)=1
do i=0,l_orbital-1
simm_orb_d(i,i,j)=1
enddo
do i=-3,3
simm_orb(i,i,j)=-1
enddo
!2: TRS TRS up=down TRS down=-up
j=2
simm_spin(0,1,j)=+1
simm_spin(1,0,j)=-1
do i=0,l_orbital-1	!d go into themselves
simm_orb_d(i,i,j)=1
enddo
do i=-3,3		!m goes into -m with a (-1)**m factor
simm_orb(i,-i,j)=(-1)**i
enddo
!M_x/i (mirror reflection)	!-i sigma_x
j=3
simm_spin(0,1,j)=-1
simm_spin(1,0,j)=-1
do i=0,l_orbital-1
simm_orb_d(i,i,j)=1
enddo
do i=-3,3
simm_orb(i,-i,j)=1
enddo
!M_y/i	!-i sigma_y 
j=4
simm_spin(0,1,j)=-(0,1)	!dosn s up=
simm_spin(1,0,j)=(0,1)
do i=0,l_orbital-1	!d
simm_orb_d(i,i,j)=1
enddo
do i=-3,3	!f
simm_orb(-i,i,j)=(-1)**i
enddo
!M_z z->-z 		!exp(-i*pi*sigma_z)=-i sigma_z
j=5
simm_spin(0,0,j)=1
simm_spin(1,1,j)=-1
do i=0,l_orbital-1	!d
simm_orb_d(i,i,j)=1
enddo
do i=-3,3	!f
simm_orb(i,i,j)=-(-1)**i
enddo
!M_xz-yz reflection along x-y plane (y=x)	!-i (sigma_x-sigma_y)
j=6
simm_spin(0,1,j)=-(1,-1)/sqrt(2.0d0)	!
simm_spin(1,0,j)=-(1,+1)/sqrt(2.0d0)
if (l_orbital>1) simm_orb_d(0,0,j)=-1	!x^2-y^2
if (l_orbital>1) simm_orb_d(1,1,j)=+1	!z^2
if (l_orbital==1)simm_orb_d(0,0,j)=+1	!s
do i=-3,3	!f
simm_orb(i,-i,j)=(0,1)**(i)
enddo
!M_xz+yz reflection along x+y plane (y=-x)	!-i (sigma_x+sigma_y) + from n=n_sf x n_mp is along -k_x,ky
j=7
simm_spin(0,1,j)=-(+1,+1)/sqrt(2.0d0)
simm_spin(1,0,j)=-(+1,-1)/sqrt(2.0d0)
if (l_orbital>1) simm_orb_d(0,0,j)=-1	!x^2-y^2
if (l_orbital>1) simm_orb_d(1,1,j)=+1	!z^2
if (l_orbital==1)simm_orb_d(0,0,j)=+1	!s
do i=-3,3	!f
simm_orb(-i,i,j)=(0,1)**(i)
enddo
!M_xz+yz reflection along x+y plane
!j=8,9,10,11
!simm_spin(0,1,j)=(-1,1)/sqrt(2.0d0)
!simm_spin(1,0,j)=(+1,1)/sqrt(2.0d0)
!do i=0,l_orbital-1	!d
!simm_orb_d(0,0,j)=-1	!x^2-y^2
!simm_orb_d(1,1,j)=+1	!z^2
!enddo
!do i=-3,3	!f
!simm_orb(-i,i,j)=(0,1)**(i)
!enddo

!R_pi/2,z 
j=8
simm_spin(0,0,j)=(1,1)/sqrt(2.0d0)	!exp (i pi/4) !exp(-pi/4 sigma_z)
simm_spin(1,1,j)=(+1,-1)/sqrt(2.0d0)	!exp (-i pi/4)
do i=0,l_orbital-1	!d
if (l_orbital>1) simm_orb_d(0,0,j)=-1	!x^2-y^2
if (l_orbital>1) simm_orb_d(1,1,j)=+1	!z^2
if (l_orbital==1)simm_orb_d(0,0,j)=+1	!z^2
enddo
do i=-3,3	!f
simm_orb(i,i,j)=(0,-1)**(i)
enddo
!R_pi,z 
j=9
simm_spin(0,0,j)=(0,1)		!exp (i pi/2)	down	!exp(-pi/2 sigma_z)
simm_spin(1,1,j)=(0,-1)		!exp (-i pi/2)  up
do i=0,l_orbital-1	!d
simm_orb_d(i,i,j)=+1	!x^2-y^2, z^2
enddo
do i=-3,3	!f
simm_orb(i,i,j)=(-1)**i
enddo
!R_pi,x 
j=10
simm_spin(0,1,j)=(0,-1)		!exp(-pi/2 sigma_x)=-i sigma_x
simm_spin(1,0,j)=(0,-1)		!
do i=0,l_orbital-1	!d
simm_orb_d(i,i,j)=+1	!x^2-y^2, z^2
enddo
do i=-3,3	!f
simm_orb(i,-i,j)=(-1)
enddo
!R_pi,x-y 
j=11
!simm_spin(0,1,j)=(1,1)/sqrt(2.0d0)	!with -pi
!simm_spin(1,0,j)=(-1,1)/sqrt(2.0d0)		
simm_spin(0,1,j)=(-1,-1)/sqrt(2.0d0)	!rotation around kx=-ky
simm_spin(1,0,j)=(1,-1)/sqrt(2.0d0)		
!simm_spin(0,1,j)=(1,-1)/sqrt(2.0d0)	!rotation around kx=ky
!simm_spin(1,0,j)=(-1,-1)/sqrt(2.0d0)		
do i=0,l_orbital-1	!d
simm_orb_d(i,i,j)=+1	!x^2-y^2
enddo
do i=-3,3	!f
simm_orb(-i,i,j)=-(0,-1)**i		!theta->pi-theta, phi->-pi/2-phi
simm_orb(-i,i,j)=-(0,1)**i		!theta->pi-theta, phi->+pi/2-phi
enddo
!M_yx-zx reflection along y-z plane (y=z)	!-i (sigma_y-sigma_z)
j=12
simm_spin(0,0,j)=(0,1)/sqrt(2.0d0)	!
simm_spin(1,1,j)=(0,-1)/sqrt(2.0d0)
simm_spin(0,1,j)=-(-1,0)/sqrt(2.0d0)	!
simm_spin(1,0,j)=-(1,0)/sqrt(2.0d0)
if (l_orbital>1) simm_orb_d(0,0,j)=0.5	!x^2-y^2
if (l_orbital>1) simm_orb_d(0,1,j)=-sqrt(3.0)/2	!z^2
if (l_orbital>1) simm_orb_d(1,0,j)=-sqrt(3.0)/2	!z^2
if (l_orbital>1) simm_orb_d(1,1,j)=-0.5	!z^2
if (l_orbital==1)simm_orb_d(0,0,j)=+1	!s
!do i=-3,3	!f
!simm_orb(i,-i,j)=(0,1)**(i)
!enddo



simm_refl=0
!f
do o1=0, l_orbital_f-1
do o2=0, l_orbital_f-1
do is=0,1	!pseudospin	
do is2=0,1	!pseudospin
do ixyz=0,n_sym_arpes
do icf=1,1
do k=-3,3	!m, to be traced
do k2=-3,3	!m, to be traced
do is3=0,1	!spin
do is4=0,1	!spin

simm_refl(is,is2,o1,o2,icf,ixyz)=simm_refl(is,is2,o1,o2,icf,ixyz)+alpha_msigma(o1,is,is3,k)*simm_spin(is3,is4,ixyz)*simm_orb(k,k2,ixyz)*alpha_msigma(o2,is2,is4,k2)
!allocate(alpha_msigma(0:l_orbital_f-1,0:1,0:1,-3:3))
!			 o1	       is  is3  k


enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo


!simm_refl(:,:,:,:,1,12)=0
!simm_refl(0,0,0,0,1,12)=(0,1)/(2*sqrt(2.0))
!simm_refl(0,1,0,0,1,12)=(-1,0)/(2*sqrt(2.0))
!if (l_orbital_f>1) simm_refl(0,0,0,1,1,12)=(0,-1)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(0,1,0,1,1,12)=(1,0)*sqrt(3./8)
!simm_refl(1,0,0,0,1,12)=(1,0)/(2*sqrt(2.0))
!simm_refl(1,1,0,0,1,12)=(0,-1)/(2*sqrt(2.0))
!if (l_orbital_f>1) simm_refl(1,0,0,1,1,12)=(-1,0)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(1,1,0,1,1,12)=(0,1)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(0,0,1,0,1,12)=(0,-1)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(0,1,1,0,1,12)=(1,0)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(0,0,1,1,1,12)=(0,-1)/(2*sqrt(2.0))
!if (l_orbital_f>1) simm_refl(0,1,1,1,1,12)=(1,0)/(2*sqrt(2.0))
!if (l_orbital_f>1) simm_refl(1,0,1,0,1,12)=(-1,0)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(1,1,1,0,1,12)=(0,1)*sqrt(3./8)
!if (l_orbital_f>1) simm_refl(1,0,1,1,1,12)=(-1,0)/(2*sqrt(2.0))
!if (l_orbital_f>1) simm_refl(1,1,1,1,1,12)=(0,1)/(2*sqrt(2.0))
!if (l_orbital_f>2) simm_refl(0,0,2,2,1,12)=(0,-1)/sqrt(2.0)
!if (l_orbital_f>2) simm_refl(0,1,2,2,1,12)=(1,0)/sqrt(2.0)
!if (l_orbital_f>2) simm_refl(0,0,2,2,1,12)=(-1,0)/sqrt(2.0)
!if (l_orbital_f>2) simm_refl(0,1,2,2,1,12)=(0,1)/sqrt(2.0)
!simm_refl(:,:,:,:,1,12)=-simm_refl(:,:,:,:,1,12)




!d
do o1=0, l_orbital-1
do o2=0, l_orbital-1
do is=0,1	!spin	
do is2=0,1	!spin
do ixyz=0,n_sym_arpes
do icf=0,0

simm_refl(is,is2,o1,o2,icf,ixyz)=simm_spin(is,is2,ixyz)*simm_orb_d(o1,o2,ixyz)

simm_refl(is,is2,o1,o2,1,12)=-simm_refl(is,is2,o1,o2,0,12)

enddo
enddo
enddo
enddo
enddo
enddo

!correction to have Hermitian operator for reflections
!do ixyz=2,7
!simm_refl(:,:,:,:,:,ixyz)=simm_refl(:,:,:,:,:,ixyz)/(0,1)
!enddo

!total
simm_oper=0
do is=0,1	!pseudospin	
do is2=0,1	!pseudospin
do ixyz=0,n_sym_arpes
do icf=0,1
do o1=0, l_orb(icf)-1
do o2=0, l_orb(icf)-1

ind=index_cso(icf,is,o1)
ind2=index_cso(icf,is2,o2)

simm_oper(ind,ind2,ixyz)=simm_refl(is,is2,o1,o2,icf,ixyz)

enddo
enddo
enddo
enddo
enddo
enddo

!squared of symm operator
simm_oper2=0
do ind=0,dim_hk-1
do ind2=0,dim_hk-1
do ind3=0,dim_hk-1
do ixyz=0,n_sym_arpes

simm_oper2(ind,ind2,ixyz)=simm_oper2(ind,ind2,ixyz)+simm_oper(ind,ind3,ixyz)*simm_oper(ind3,ind2,ixyz)

enddo
enddo
enddo
enddo


!M^2=-1 checks for consistency
do ind=0,dim_hk-1
do ind2=0,dim_hk-1
do ixyz=0,n_sym_arpes

if (ind==ind2) then
if (abs(abs(simm_oper2(ind,ind2,ixyz))-1)>1e-5) print *, "WARNING refl_symm",  ixyz
else
if (abs(abs(simm_oper2(ind,ind2,ixyz)))>1e-5) print *, "WARNING refl_symm",  ixyz
endif

enddo
enddo
enddo



print *, "symmetries"


do ixyz=0, n_sym_arpes
print *, "Sym n.:", ixyz

do i=0, dim_hk-1
!write (*, "(10(f.1,f5.1))"), (simm_oper(i,j,ixyz), j=0, dim_hk-1)
write(*,"(10(a1,f5.2,a1,f5.2,a1))") ( "(",real(simm_oper(i,j,ixyz)),",",aimag(simm_oper(i,j,ixyz)),")", j=0,dim_hk-1 )
enddo

print *, ""

enddo

!do ixyz=0, n_sym_arpes
!print *, "Sym n.:", ixyz

!do i=0, dim_hk-1
!write (*, "(10(f.1,f5.1))"), (simm_oper(i,j,ixyz), j=0, dim_hk-1)
!write(*,"(10(a1,f5.2,a1,f5.2,a1))") ( "(",real(simm_oper2(i,j,ixyz)),",",aimag(simm_oper2(i,j,ixyz)),")", j=0,dim_hk-1 )
!enddo

!print *, ""

!enddo




!simm_refl=0

!do o1=0, l_orbital_f-1
!Identity
!simm_refl(0,0,o1,o1,:,0)=1
!simm_refl(1,1,o1,o1,:,0)=1
!M_yz
!simm_refl(0,1,o1,o1,:,1)=(0,1)
!simm_refl(1,0,o1,o1,:,1)=(0,1)
!M_xz
!simm_refl(0,1,o1,o1,0,2)=-1		!d
!simm_refl(1,0,o1,o1,0,2)=+1
!simm_refl(0,1,o1,o1,1,2)=+1		!f
!simm_refl(1,0,o1,o1,1,2)=-1
!M_xy
!simm_refl(0,0,o1,o1,0,3)=-(0,1)		!d
!simm_refl(1,1,o1,o1,0,3)=+(0,1)
!simm_refl(0,0,o1,o1,1,3)=+(0,1)		!f
!simm_refl(1,1,o1,o1,1,3)=-(0,1)
!M_(x-y)z
!simm_refl(0,1,o1,o1,0,4)=(1,1)/sqrt(2.)		!d
!simm_refl(1,0,o1,o1,0,4)=(-1,1)/sqrt(2.)
!if (o1==0) simm_refl(0,1,o1,o1,1,4)=(1,1)/sqrt(2.)
!if (o1==0) simm_refl(1,0,o1,o1,1,4)=(-1,1)/sqrt(2.)
!if (o1==1) simm_refl(0,1,o1,o1,1,4)=(1,1)/sqrt(2.)
!if (o1==1) simm_refl(1,0,o1,o1,1,4)=(-1,1)/sqrt(2.)
!time reversal
!simm_refl(0,0,o1,o1,0,3)=-(0,1)		!d
!simm_refl(1,1,o1,o1,0,3)=+(0,1)
!simm_refl(0,0,o1,o1,1,3)=+(0,1)		!f
!simm_refl(1,1,o1,o1,1,3)=-(0,1)
!enddo

!sigma_refl(is,is2,o1,o2,icf,ixyz)


!print *, "Symmetries on Basis:"
!print *, "s-d"
!icf=0

!if (l_orbital==2) then

!ixyz=1
!write (*, "(a5, 4a18)"), "M_yz", "d1-", "d1+", "d2-", "d2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)

!ixyz=2
!write (*, "(a5, 4a18)"), "M_xz", "d1-", "d1+", "d2-", "d2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)

!ixyz=3
!write (*, "(a5, 4a18)"), "M_xy", "d1-", "d1+", "d2-", "d2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)

!print *, ""

!else

!ixyz=1
!write (*, "(a5, 4a18)"), "M_yz", "d-", "d+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz)

!ixyz=2
!write (*, "(a5, 4a18)"), "M_xz", "d-", "d+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz)

!ixyz=3
!write (*, "(a5, 4a18)"), "M_xy", "d-", "d+"
!write (*, "(a5, 4(f10.4,f8.4))"), "d-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "d+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz)

!print *, ""

!endif


!print *, "f"
!icf=1

!if (l_orbital_f==2) then

!ixyz=1
!write (*, "(a5, 4a18)"), "M_yz", "f1-", "f1+", "f2-", "f2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)

!ixyz=2
!write (*, "(a5, 4a18)"), "M_xz", "d1-", "d1+", "d2-", "d2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)

!ixyz=3
!write (*, "(a5, 4a18)"), "M_xy", "f1-", "f1+", "f2-", "f2+"
!write (*, "(a5, 4(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz)
!write (*, "(a5, 4(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz)


!else if (l_orbital_f==3) then
!ixyz=1
!write (*, "(a5, 6a18)"), "M_yz", "f1-", "f1+", "f2-", "f2+", "f3-", "f3+"
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz),&
!					 simm_refl(0,0,0,2,icf,ixyz),simm_refl(0,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz),&
!					 simm_refl(1,0,0,2,icf,ixyz),simm_refl(1,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz),&
!					 simm_refl(0,0,1,2,icf,ixyz),simm_refl(0,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz),&
!					 simm_refl(1,0,1,2,icf,ixyz),simm_refl(1,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,2,0,icf,ixyz),simm_refl(0,1,2,0,icf,ixyz),simm_refl(0,0,2,1,icf,ixyz),simm_refl(0,1,2,1,icf,ixyz),&
!					 simm_refl(0,0,2,2,icf,ixyz),simm_refl(0,1,2,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,2,0,icf,ixyz),simm_refl(1,1,2,0,icf,ixyz),simm_refl(1,0,2,1,icf,ixyz),simm_refl(1,1,2,1,icf,ixyz),&
!					 simm_refl(1,0,2,2,icf,ixyz),simm_refl(1,1,2,2,icf,ixyz)

!ixyz=2
!write (*, "(a5, 4a18)"), "M_xz", "f1-", "f1+", "f2-", "f2+","f3-", "f3+"
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz),&
!					 simm_refl(0,0,0,2,icf,ixyz),simm_refl(0,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz),&
!					 simm_refl(1,0,0,2,icf,ixyz),simm_refl(1,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz),&
!					 simm_refl(0,0,1,2,icf,ixyz),simm_refl(0,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz),&
!					 simm_refl(1,0,1,2,icf,ixyz),simm_refl(1,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,2,0,icf,ixyz),simm_refl(0,1,2,0,icf,ixyz),simm_refl(0,0,2,1,icf,ixyz),simm_refl(0,1,2,1,icf,ixyz),&
!					 simm_refl(0,0,2,2,icf,ixyz),simm_refl(0,1,2,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,2,0,icf,ixyz),simm_refl(1,1,2,0,icf,ixyz),simm_refl(1,0,2,1,icf,ixyz),simm_refl(1,1,2,1,icf,ixyz),&
!					 simm_refl(1,0,2,2,icf,ixyz),simm_refl(1,1,2,2,icf,ixyz)
!ixyz=3
!write (*, "(a5, 4a18)"), "M_xy", "f1-", "f1+", "f2-", "f2+", "f3-", "f3+"
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,0,0,icf,ixyz),simm_refl(0,1,0,0,icf,ixyz),simm_refl(0,0,0,1,icf,ixyz),simm_refl(0,1,0,1,icf,ixyz),&
!					 simm_refl(0,0,0,2,icf,ixyz),simm_refl(0,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,0,0,icf,ixyz),simm_refl(1,1,0,0,icf,ixyz),simm_refl(1,0,0,1,icf,ixyz),simm_refl(1,1,0,1,icf,ixyz),&
!					 simm_refl(1,0,0,2,icf,ixyz),simm_refl(1,1,0,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2-", simm_refl(0,0,1,0,icf,ixyz),simm_refl(0,1,1,0,icf,ixyz),simm_refl(0,0,1,1,icf,ixyz),simm_refl(0,1,1,1,icf,ixyz),&
!					 simm_refl(0,0,1,2,icf,ixyz),simm_refl(0,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f2+", simm_refl(1,0,1,0,icf,ixyz),simm_refl(1,1,1,0,icf,ixyz),simm_refl(1,0,1,1,icf,ixyz),simm_refl(1,1,1,1,icf,ixyz),&
!					 simm_refl(1,0,1,2,icf,ixyz),simm_refl(1,1,1,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1-", simm_refl(0,0,2,0,icf,ixyz),simm_refl(0,1,2,0,icf,ixyz),simm_refl(0,0,2,1,icf,ixyz),simm_refl(0,1,2,1,icf,ixyz),&
!					 simm_refl(0,0,2,2,icf,ixyz),simm_refl(0,1,2,2,icf,ixyz)
!write (*, "(a5, 6(f10.4,f8.4))"), "f1+", simm_refl(1,0,2,0,icf,ixyz),simm_refl(1,1,2,0,icf,ixyz),simm_refl(1,0,2,1,icf,ixyz),simm_refl(1,1,2,1,icf,ixyz),&
!					 simm_refl(1,0,2,2,icf,ixyz),simm_refl(1,1,2,2,icf,ixyz)

!endif

print *, "*******"







end subroutine set_spin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_xyzcso(ix,iy,iz,icf,is,o1,nx,ny,nz)
integer, intent(in):: ix,iy,iz,icf,is,o1,nz,nx,ny

if (l_orbital==l_orbital_f .or. o1 < l_orbital) then
 index_xyzcso =ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz + o1*4*nx*ny*nz
!elseif (l_orbital==0 .and. l_orbital_f==1) then
!index_xyzcso =ix + iy*nx + iz*nx*ny + icf*nx*ny*nz + is*2*nx*ny*nz + o1*4*nx*ny*nz
else
 index_xyzcso =ix + iy*nx + iz*nx*ny + (icf-1)*nx*ny*nz + is*nx*ny*nz + o1*4*nx*ny*nz
endif
if (l2dmodel) index_xyzcso=ix + iy*nx + iz*nx*ny + is*nx*ny*nz

if (index_xyzcso<0) write (*, "(a13,6i4, i10)"), "index_xyzcso",ix,iy,iz,icf,is,o1,index_xyzcso
!write (*, "(i10)"), index_xyzcso

end function index_xyzcso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_xyzcso2(ix,iy,iz,j,nx,ny,nz)
integer, intent(in):: ix,iy,iz,j,nz,nx,ny

 index_xyzcso2 =ix + iy*nx + iz*nx*ny + j*nx*ny*nz
end function index_xyzcso2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_xyz(ix,iy,iz,nx,ny)
integer, intent(in):: ix,iy,iz,nx,ny

index_xyz=ix+nx*iy+nx*ny*iz

end function index_xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_cso(icf,is,o1)
integer, intent(in):: icf,is,o1

!c1-c1+c2-c2+f1-f1+f2-f2+f3-f3+

 index_cso = is+2*o1+2*(l_orbital)*icf
if (l2dmodel) index_cso=is

end function index_cso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_zcso(iz,icf,is,o1,nz)
integer, intent(in):: iz,icf,is,o1,nz

if (l_orbital==l_orbital_f .or. o1 < l_orbital) then
 index_zcso =iz + icf*nz + is*2*nz + o1*4*nz
elseif (l_orbital_f==1 .and. l_orbital==0) then !for l2dmodel
index_zcso=is
else
 index_zcso =iz + (icf-1)*nz + is*nz + o1*4*nz
endif

if (l2dmodel) index_zcso=is

end function index_zcso


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_xyzcso_red(ix,iy,ix2,iy2,iz,icf,is,o1,nz, nx_red, nxy_red)
integer, intent(in):: ix,iy,iz,icf,is,o1,nz, nx_red, nxy_red,ix2,iy2

!ix,iy= index of cell
!ix2,iy2= index of sites inside cell

if(lc4v) then	!Here I assume nx_kr=1 and ny_kr=1, so don't use ix2,iy2 which are zero

if (l_orbital==l_orbital_f .or. o1 < l_orbital) then
 index_xyzcso_red = ix + iy*nx_red - (iy*(iy+1))/2 + iz*nxy_red + icf*nxy_red*nz_plot + is*2*nxy_red*nz_plot + o1*4*nxy_red*nz_plot
else
 index_xyzcso_red = ix + iy*nx_red - (iy*(iy+1))/2 + iz*nxy_red + (icf-1)*nxy_red*nz_plot + is*nxy_red*nz_plot + o1*4*nxy_red*nz_plot
endif

else if(lc2v .or. lcs .or. lnos) then!c2v symmetry I think I assume nx=nx_kr but I am not sure

if (l_orbital==l_orbital_f .or. o1 < l_orbital) then
 index_xyzcso_red = ix2 + ix*nx_kr + iy2*nx_red*nx_kr + iy*ny_kr*nx_red*nx_kr + iz*nxy_red*nx_kr*ny_kr + icf*nxy_red*nz_plot*nx_kr*ny_kr + &
 		    is*2*nxy_red*nz_plot*nx_kr*ny_kr + o1*4*nxy_red*nz_plot*nx_kr*ny_kr
else
 index_xyzcso_red = ix2 + ix*nx_kr + iy2*nx_red*nx_kr + iy*ny_kr*nx_red*nx_kr + iz*nxy_red*nx_kr*ny_kr + (icf-1)*nxy_red*nz_plot*nx_kr*ny_kr + &
 		    is*nxy_red*nz_plot*nx_kr*ny_kr + o1*4*nxy_red*nz_plot*nx_kr*ny_kr
endif

!else if(lcs) then !c_s symmetry I think I assume nx=nx_kr but I am not sure

!if (l_orbital==l_orbital_f .or. o1 < l_orbital) then
! index_xyzcso_red = ix2 + ix*nx_kr + iy2*nkx_plot*nx_kr + iy*ny_kr*nkx_plot*nx_kr + iz*nxy_red*nx_kr*ny_kr + icf*nxy_red*nz_plot*nx_kr*ny_kr + &
! 		    is*2*nxy_red*nz_plot*nx_kr*ny_kr + o1*4*nxy_red*nz_plot*nx_kr*ny_kr
!else
! index_xyzcso_red = ix2 + ix*nx_kr + iy2*nkx_plot*nx_kr + iy*ny_kr*nkx_plot*nx_kr + iz*nxy_red*nx_kr*ny_kr + (icf-1)*nxy_red*nz_plot*nx_kr*ny_kr + &
 !		    is*nxy_red*nz_plot*nx_kr*ny_kr + o1*4*nxy_red*nz_plot*nx_kr*ny_kr	!nxy_red=nkx_plot*(nky_plot/2+1), nkx_plot=nx_red, ny_red=nky_plot/2+1
!endif



else
print *, "wronf lc4v, lc2v, lcs"
stop

endif

if (l2dmodel) index_xyzcso_red=ix + iy*nx_red - (iy*(iy+1))/2 + is*nxy_red*nz_plot

!write (*, "(6i4, i10)"), ix,iy,iz,icf,is,o1,index_xyzcso
!write (*, "(i10)"), index_xyzcso

end function index_xyzcso_red



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_xyz_red(ix,iy,ix2,iy2,iz,nz, nx_red, nxy_red)
integer, intent(in):: ix,iy,iz,nz, nx_red, nxy_red,ix2,iy2

!ix=index of cell
!ix2=index of atom in cell

if(lc4v) then	!Here I assume nx_kr=1 and ny_kr=1, so don't use ix2,iy2 which are zero
index_xyz_red = ix + iy*nx_red - (iy*(iy+1))/2 + iz*nxy_red
elseif(lc2v .or. lcs .or. lnos) then !c2v symmetry
 index_xyz_red = ix2 + ix*nx_kr + iy2*nx_red*nx_kr + iy*ny_kr*nx_red*nx_kr + iz*nxy_red*nx_kr*ny_kr
!elseif(lcs) then
! index_xyz_red = ix2 + ix*nx_kr + iy2*nkx_plot*nx_kr + iy*ny_kr*nx_kr*nkx_plot + iz*nxy_red*nx_kr*ny_kr
endif

if (l2dmodel) index_xyz_red=ix + iy*nx_red - (iy*(iy+1))/2


end function index_xyz_red

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function index_co(icf,o1)
integer, intent(in):: icf,o1

!c1c2f1f2f3
!c1c2f1f2
!c1f1

 index_co = o1+l_orbital*icf

end function index_co
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readham 

allocate(Ham_input_temp(-n_shells:n_shells, -n_shells:n_shells,-n_shells:n_shells, 0:dim_hk-1, 0:dim_hk-1))
allocate(Ham_input(-n_shells:n_shells, -n_shells:n_shells,-n_shells:n_shells, 0:dim_hk-1, 0:dim_hk-1))
Ham_input_temp=0
Ham_input=0

OPEN(unit=101,file='MLWF_hr.dat',status='unknown')

print *, "Reading Ham"
print *, "n_shells", n_shells
!print *, "n_entry", n_entry_in
print *, "dist", dist
print *, "dist_max", dist_max


do l=1,n_entry_in
 read (101, format_in) i,j,k,m,n, real_h, imag_h


if ( (abs(i)< n_shells+1) .and. (abs(j)< n_shells+1) .and. (abs(k)< n_shells+1) .and. (m < dim_hk+1) .and. (n < dim_hk+1)) &
Ham_input_temp(i,j,k,m-1,n-1)=-real_h-(0,1.0)*imag_h		!I take minus because I use holes

enddo


 close(unit=101)

print *, "Reading OK"

print *, "Rotating Ham"
!into gamma 78 basis and x2y2 down x2y2 up z2 down z2 up 

allocate(Ug78(0:dim_hk-1,0:dim_hk-1))
allocate(Ug78m1(0:dim_hk-1,0:dim_hk-1))

call set_Ug78 (Ug78, Ug78m1)

do ix=-n_shells, n_shells
do iy=-n_shells, n_shells
do iz=-n_shells, n_shells

do m=0, dim_hk-1
do n=0, dim_hk-1
do i=0, dim_hk-1
do j=0, dim_hk-1

Ham_input(ix,iy,iz,m,n)=Ham_input(ix,iy,iz,m,n)+Ug78(m,i)*Ham_input_temp(ix,iy,iz,i,j)*Ug78m1(j,n)


enddo
enddo
enddo
enddo
enddo
enddo
enddo


!shift f energies
do i=4,7
Ham_input(0,0,0,i,i)=Ham_input(0,0,0,i,i)+delta_e8_read
enddo

do i=8,9
Ham_input(0,0,0,i,i)=Ham_input(0,0,0,i,i)+delta_e7_read
enddo

!rescales f kin energy
do ix=-n_shells, n_shells
do iy=-n_shells, n_shells
do iz=-n_shells, n_shells

do i=4,9
do j=4,9

if (ix .ne. 0 .or. iy .ne. 0 .or. iz .ne. 0) Ham_input(ix,iy,iz,i,j)=Ham_input(ix,iy,iz,i,j)*delta_tf_read


enddo
enddo
enddo
enddo
enddo

!rescales c kin energy
do ix=-n_shells, n_shells
do iy=-n_shells, n_shells
do iz=-n_shells, n_shells

do i=0,3
do j=0,3

if (ix .ne. 0 .or. iy .ne. 0 .or. iz .ne. 0) Ham_input(ix,iy,iz,i,j)=Ham_input(ix,iy,iz,i,j)*delta_tc_read


enddo
enddo
enddo
enddo
enddo


!rescales hybridization
do ix=-n_shells, n_shells
do iy=-n_shells, n_shells
do iz=-n_shells, n_shells

do i=4,9
do j=0,3

if (ix .ne. 0 .or. iy .ne. 0 .or. iz .ne. 0) Ham_input(ix,iy,iz,i,j)=Ham_input(ix,iy,iz,i,j)*delta_v_read
if (ix .ne. 0 .or. iy .ne. 0 .or. iz .ne. 0) Ham_input(ix,iy,iz,j,i)=Ham_input(ix,iy,iz,j,i)*delta_v_read


enddo
enddo
enddo
enddo
enddo


deallocate(Ham_input_temp)
deallocate(Ug78)
deallocate(Ug78m1)

print *, "Rotating OK"

 
 
 
end subroutine readham 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_Ug78(U, Um1)
complex(kind=idpc), intent(out)	::U(0:dim_hk-1,0:dim_hk-1),Um1(0:dim_hk-1,0:dim_hk-1)
complex(kind=idpc)	::UstarU(0:dim_hk-1,0:dim_hk-1)
integer		::ind,ix,iy,iz,ik, i,j,k

print *, "Setting U"

U=0

!d: x2-y2 down, x2-y2 up, z2 down, z2 up
U(0,3)=1
U(1,2)=1
U(2,1)=1
U(3,0)=1

!Gamma81-
U(4,4)=sqrt(5./6.)
U(4,8)=sqrt(1./6.)
!Gamma81+
U(5,5)=sqrt(1./6.)
U(5,9)=sqrt(5./6.)
!Gamma82-
U(6,6)=1.
!Gamma82+
U(7,7)=1.
!Gamma7-
U(8,4)=+sqrt(1./6.)
U(8,8)=-sqrt(5./6.)
!Gamma7+
U(9,5)=-sqrt(5./6.)
U(9,9)=+sqrt(1./6.)


!f: gamma81-, gamma81+, gamma82-, gamma82+, gamma7-, gamma7+



do i=0,dim_hk-1
  do j=0,dim_hk-1

   Um1(i,j)=conjg(U(j,i))
 
  enddo
enddo


UstarU=0

do i=0,dim_hk-1
  do j=0,dim_hk-1
     do k=0,dim_hk-1

   UstarU (i,j)=UstarU(i,j)+Um1(i,k)*U(k,j)
   
     enddo
  enddo
enddo
	
	
!call write_ham_k(UstarU)

!call write_ham(UstarU,"UstarU",dim_hk,dim_hk)


do i=0,dim_hk-1
  do j=0,dim_hk-1

   if (i==j .and. abs(Ustaru(i,j)-1)>1e-6) then
   print *, "error in U",i, j , Ustaru(i,j)
 !  stop
   endif

   if (i.ne.j .and. abs(Ustaru(i,j))>1e-6) then
   print *, "error in U", i, j , Ustaru(i,j)
!   stop
   endif

  enddo
enddo



print *, "OK setting U"
	
end subroutine set_Ug78
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp) function mod2(a,b)
real(kind=idp):: a,b

if(a<=b/2) mod2=a
if(a>b/2) mod2=a-b

return

end function mod2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_ixmax(ixmax)
integer,intent(out)::ixmax

ixmax=nx_red-1
!if(lc4v) ixmax=nx_red-1
!if(lc2v) ixmax=nx_red-1
!if(lcs)  ixmax=nkx_plot-1	!but here nkx_plot=nx_red

end subroutine set_ixmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_iymax(iymax,ix)
integer,intent(out)::iymax
integer,intent(in)::ix

if(lc4v) iymax=ix
if(lc2v .or. lcs .or. lnos) iymax=ny_red-1
!if(lcs)  iymax=ny_red-1

end subroutine set_iymax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine index_qk1k2(iq,ik1,ik2,nkx,nky)
integer,intent(out)::iq
integer,intent(in)::ik1,ik2,nkx,nky
integer:: ik1x,ik1y,ik2x,ik2y,iqx,iqy, ikx0,iky0


!iqx=mod(ikx-ikx2+3*nkx_plot/2,nkx_plot)
!iqy=mod(iky-iky2+3*nkx_plot/2,nkx_plot)


if(mod(nkx,2)==0) ikx0=nkx/2
if(mod(nkx,2)==1) ikx0=(nkx-1)/2
if(mod(nky,2)==0) iky0=nky/2
if(mod(nky,2)==1) iky0=(nky-1)/2

ik1x=mod(ik1,nkx)
ik1y=(ik1-ik1x)/nkx
ik2x=mod(ik2,nkx)
ik2y=(ik2-ik2x)/nkx
iqx=mod(ik1x-ik2x+3*nkx/2,nkx)
iqy=mod(ik1y-ik2y+3*nky/2,nky)

iq=iqx+iqy*nkx

!write(*, "(10i6)"),ik1,ik1x,ik1y,ik2,ik2x,ik2y,iq,iqx,iqy


end subroutine index_qk1k2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ik_pxpy (ik,nkx, nky, ik_px, ik_py)	!given ik, it gives ik_px just on the right, and ik_py just above (within pbc)
integer,intent(in)::ik,nkx,nky
integer,intent(out)::ik_px,ik_py
integer::ikx,iky,ikx_px,iky_px,ikx_py,iky_py


ikx=mod(ik,nkx)		!ik=ikx+nkx*iky
iky=(ik-ikx)/nkx

!print *,"i", ik,ikx,iky

iky_px=iky
if(ikx==nkx-1) then
ikx_px=0
else
ikx_px=ikx+1
endif

ikx_py=ikx
if(iky==nky-1) then
iky_py=0
else
iky_py=iky+1
endif


ik_px=ikx_px+nkx*iky_px
ik_py=ikx_py+nkx*iky_py

!print *, ik, ik_px, ik_py



end subroutine ik_pxpy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind=idpc) function det(N, mat)

integer, intent(in) :: N 
complex(kind=idpc), intent(inout), dimension(0:n-1,0:n-1) :: mat
integer:: i, info
integer, allocatable :: ipiv(:)
real(kind=idp) :: sgn

allocate(ipiv(0:n-1))
ipiv = 0

call zgetrf(N, N, mat, N, ipiv, info)

det = 1.0d0

do i = 0, N-1
  det = det*mat(i, i)
end do

sgn = 1.0d0

do i = 0, N-1
 if(ipiv(i) /= i) then
  sgn = -sgn
 end if
 end do

 det = sgn*det   


deallocate(ipiv)

!if (n==2) det=mat(0,0)*mat(1,1)-mat(0,1)*mat(1,0)
!if (n==3) det=mat(0,0)*(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))-mat(0,1)*(mat(1,0)*mat(2,2)-mat(1,2)*mat(2,0))+mat(0,2)*(mat(1,0)*mat(2,1)-mat(1,1)*mat(2,0))
!I need the determinant :-(

!call diagonalize_s(n_occ/2,u_mat_fm,e_u_fm)
!call diagonalize_s(n_occ/2,u_mat_fp,e_u_fp)
!call diagonalize_s(n_empty/2,u_mat_em,e_u_em)
!call diagonalize_s(n_empty/2,u_mat_ep,e_u_ep)

!u_k_chern(ik,1,1)=product(e_u_fm)/abs(product(e_u_fm))
!u_k_chern(ik,1,2)=product(e_u_fp)/abs(product(e_u_fp))
!u_k_chern(ik,1,3)=product(e_u_em)/abs(product(e_u_em))
!u_k_chern(ik,1,4)=product(e_u_ep)/abs(product(e_u_ep))



end function det
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp) function det_real(N, mat)

integer, intent(in) :: N 
real(kind=idp), intent(inout), dimension(0:n-1,0:n-1) :: mat
integer:: i, info
integer, allocatable :: ipiv(:)
real(kind=idp) :: sgn

allocate(ipiv(0:n-1))
ipiv = 0

call zgetrf(N, N, mat, N, ipiv, info)

det_real = 1.0d0

do i = 0, N-1
  det_real = det_real*mat(i, i)
end do

sgn = 1.0d0

do i = 0, N-1
 if(ipiv(i) /= i) then
  sgn = -sgn
 end if
 end do

 det_real = sgn*det_real   


deallocate(ipiv)

!if (n==2) det=mat(0,0)*mat(1,1)-mat(0,1)*mat(1,0)
!if (n==3) det=mat(0,0)*(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))-mat(0,1)*(mat(1,0)*mat(2,2)-mat(1,2)*mat(2,0))+mat(0,2)*(mat(1,0)*mat(2,1)-mat(1,1)*mat(2,0))
!I need the determinant :-(

!call diagonalize_s(n_occ/2,u_mat_fm,e_u_fm)
!call diagonalize_s(n_occ/2,u_mat_fp,e_u_fp)
!call diagonalize_s(n_empty/2,u_mat_em,e_u_em)
!call diagonalize_s(n_empty/2,u_mat_ep,e_u_ep)

!u_k_chern(ik,1,1)=product(e_u_fm)/abs(product(e_u_fm))
!u_k_chern(ik,1,2)=product(e_u_fp)/abs(product(e_u_fp))
!u_k_chern(ik,1,3)=product(e_u_em)/abs(product(e_u_em))
!u_k_chern(ik,1,4)=product(e_u_ep)/abs(product(e_u_ep))



end function det_real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind=idpc) function det_man(N, mat)
integer, intent(in) :: N 
complex(kind=idpc), intent(in), dimension(0:n-1,0:n-1) :: mat

if (n==0) det_man=0
if (n==1) det_man=mat(0,0)
if (n==2) det_man=mat(0,0)*mat(1,1)-mat(0,1)*mat(1,0)
if (n==3) det_man=mat(0,0)*(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))-mat(0,1)*(mat(1,0)*mat(2,2)-mat(1,2)*mat(2,0))+mat(0,2)*(mat(1,0)*mat(2,1)-mat(1,1)*mat(2,0))



end function det_man

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine index_01(parity,n_states,ind_nu)
integer,intent(in)::ind_nu,n_states
real(kind=idp), intent(in) :: parity(0:7, 0:n_states-1)
integer :: n_p(0:7), n_m(0:7), prod_parities, ind_01
logical	:: ok

n_p=0
n_m=0
ok=.true.

do ik=0,7

 do j=0,n_states-1
 if (abs(parity(ik,j)-1.0d0)<1e-5) then
 n_p(ik)=n_p(ik)+1
 elseif (abs(parity(ik,j)+1.0d0)<1e-5)  then
 n_m(ik)=n_m(ik)+1
 else
 print *, "error in parity product", parity(ik,j)
 endif
 enddo

!check if they are even
 if(mod(n_p(ik),2)==1 .or. mod(n_m(ik),2)==1) then
  print *, "odd number of minus or plus in ik=", ik
 ok=.false.
 endif
enddo

!write (*, "(10i5)"), (n_m(j), j=0,7)


!if they are even, divide by 2
! n_p(:)=n_p(:)/2
 n_m(:)=n_m(:)/2

!write (*, "(10i5)"), (n_m(j), j=0,7)
!write (*, "(10i5)"), (mod(n_m(j),2), j=0,7)

!transform from sum of minus states to -1 or +1
do ik=0,7
if(mod(n_m(ik),2)==0) then
n_m(ik)=1
else if(mod(n_m(ik),2)==1) then
 n_m(ik)=-1
endif
enddo

!write (*, "(10i5)"), (n_m(j), j=0,7)


!product of parities at hsp
!0(0,0,0)
!1(pi,0,0)
!2(0,pi,0)
!3(0,0,pi)
!4(pi,pi,0)
!5(pi,0,pi)
!6(0,pi,pi)
!7(pi,pi,pi)
if (ind_nu==0) prod_parities=product(n_m)
if (ind_nu==1) prod_parities=n_m(1)*n_m(4)*n_m(5)*n_m(7)
if (ind_nu==2) prod_parities=n_m(2)*n_m(4)*n_m(6)*n_m(7)
if (ind_nu==3) prod_parities=n_m(3)*n_m(5)*n_m(6)*n_m(7)
!if (ind_nu==1) prod_parities=n_m(0)*n_m(2)*n_m(3)*n_m(6)
!if (ind_nu==2) prod_parities=n_m(0)*n_m(1)*n_m(3)*n_m(5)
!if (ind_nu==3) prod_parities=n_m(0)*n_m(1)*n_m(2)*n_m(4)

!print *, prod_parities

!from -1,1 to 1,0
if (prod_parities==-1) ind_01=1
if (prod_parities==+1)  ind_01=0
if (prod_parities .ne. -1 .and. prod_parities .ne. +1) ok=.false.

if (.not. ok) ind_01=-100


write (*, "(a3,i1,i5)"), "vu_",ind_nu,ind_01
write (300, "(a3,i1,i5)"), "vu_",ind_nu,ind_01


end subroutine index_01
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function x_ss(ix,iy)
integer, intent(in)	:: ix,iy

x_ss=x_surf(1,1)*ix+x_surf(2,1)*iy
end function x_ss
real function y_ss(ix,iy)
integer, intent(in)	:: ix,iy

y_ss=x_surf(1,2)*ix+x_surf(2,2)*iy
end function y_ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function det3_real(r1r2r3)
real(kind=idp),intent(in) ::r1r2r3(3,3)

det3_real=r1r2r3(1,1)*(r1r2r3(2,2)*r1r2r3(3,3)-r1r2r3(2,3)*r1r2r3(3,2))-r1r2r3(1,2)*(r1r2r3(2,1)*r1r2r3(3,3)-r1r2r3(2,3)*r1r2r3(3,1))+r1r2r3(1,3)*(r1r2r3(2,1)*r1r2r3(3,2)-r1r2r3(2,2)*r1r2r3(3,1))
end function det3_real

integer function det3_int(r1r2r3)
integer,intent(in) ::r1r2r3(3,3)

det3_int=r1r2r3(1,1)*(r1r2r3(2,2)*r1r2r3(3,3)-r1r2r3(2,3)*r1r2r3(3,2))-r1r2r3(1,2)*(r1r2r3(2,1)*r1r2r3(3,3)-r1r2r3(2,3)*r1r2r3(3,1))+r1r2r3(1,3)*(r1r2r3(2,1)*r1r2r3(3,2)-r1r2r3(2,2)*r1r2r3(3,1))
end function det3_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_kvecs_cubic(k_vecs_nscf,nk_plot,nkx_plot)
integer,intent(in)::nkx_plot
integer,intent(out)::nk_plot
real(kind=idp),allocatable, intent(out)	::k_vecs_nscf(:,:)
integer:: m,i,j,k,l,n,ik, total_weight



m=nkx_plot/2+1	!number of independent points along x

nk_plot=0
do i=0,m
do j=0,i
do k=0,j
nk_plot=nk_plot+1
enddo
enddo
enddo

allocate(k_vecs_nscf(0:nk_plot-1,4))
k_vecs_nscf=0

ik=0
do i=0,m-1
do j=0,i
do k=0,j
if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,1)=dble(i)/dble(m-1)*0.50d0
if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,2)=dble(j)/dble(m-1)*0.50d0
if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,3)=dble(k)/dble(m-1)*0.50d0
!if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,1)=dble(i)/dble(m)*0.50d0
!if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,2)=dble(j)/dble(m)*0.50d0
!if (mod(nkx_plot,2)==0) k_vecs_nscf(ik,3)=dble(k)/dble(m)*0.50d0

if (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. i.ne.j .and. j.ne. k .and. k.ne.i) then
k_vecs_nscf(ik,4)=48	!xyz
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1)==0 .and. i.ne.j)  then
k_vecs_nscf(ik,4)=24	!xy0
elseif (mod(i,m-1)==0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. j.ne.k)  then
k_vecs_nscf(ik,4)=24	!0yz
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1)==0 .and. mod(k,m-1).ne.0 .and. i.ne.k)  then
k_vecs_nscf(ik,4)=24	!x0z
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. i==j .and. j.ne.k)  then
k_vecs_nscf(ik,4)=24	!xxz
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. j==k .and. k.ne.i)  then
k_vecs_nscf(ik,4)=24	!xyy
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. k==i .and. i.ne.j) then
k_vecs_nscf(ik,4)=24	!xyx
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. j.ne.k)  then
k_vecs_nscf(ik,4)=12	!x00.5
elseif (mod(i,m-1)==0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1)==0 .and. i.ne.k)  then
k_vecs_nscf(ik,4)=12	!0x0.5
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1).ne.0 .and. i.ne.j)  then
k_vecs_nscf(ik,4)=12	!00.5x
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1)==0 .and. i==j)  then
k_vecs_nscf(ik,4)=12	!xx0
elseif (mod(i,m-1)==0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. j==k)  then
k_vecs_nscf(ik,4)=12	!0xx
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1)==0 .and. mod(k,m-1).ne.0 .and. i==k)  then
k_vecs_nscf(ik,4)=12	!x0x
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1).ne.0 .and. i==j.and. j==k)  then
k_vecs_nscf(ik,4)=8	!xxx
elseif (mod(i,m-1).ne.0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. j==k)  then
k_vecs_nscf(ik,4)=6	!x00
elseif (mod(i,m-1)==0 .and. mod(j,m-1).ne.0 .and. mod(k,m-1)==0 .and. i==k)  then
k_vecs_nscf(ik,4)=6	!0x0
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1).ne.0 .and. i==j)  then
k_vecs_nscf(ik,4)=6	!00x
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. i==j .and. j.ne.k)  then
k_vecs_nscf(ik,4)=3	!(00 0.5)
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. i.ne.j .and. j==k)  then
k_vecs_nscf(ik,4)=3	!(0.5 00)
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. i.ne.j .and. i==k)  then
k_vecs_nscf(ik,4)=3	!(0 0.5 0)
elseif (mod(i,m-1)==0 .and. mod(j,m-1)==0 .and. mod(k,m-1)==0 .and. i==j .and. j==k) then
k_vecs_nscf(ik,4)=1	!(000) (0.5 0.5 0.5)
endif

ik=ik+1
enddo
enddo
enddo

total_weight=sum(k_vecs_nscf(:,4))

k_vecs_nscf(:,4)=k_vecs_nscf(:,4)/total_weight

print *, "total weight=", total_weight

end subroutine set_kvecs_cubic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp) function norm_k(k,nz)
integer, intent(in)::k,nz
integer::i

norm_k=0

do i=1,nz
norm_k=norm_k+(sin(pi*i*k/(nz+1)))**2
enddo

norm_k=1/sqrt(norm_k)


end function norm_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<k2|V|k>
real(kind=idp) function hybr_k(k2,k,nz)
integer, intent(in)::k,k2,nz
integer::i

hybr_k=0

if (k2 .ne. k) then
do i=1,nz
!hybr_k=hybr_k+sin(pi*i*k2/(nz+1))*cos(pi*i*k/(nz+1))
hybr_k=hybr_k+sin(pi*i*k2/(nz+1))*(sin(pi*(i-1)*k/(nz+1))-sin(pi*(i+1)*k/(nz+1)))/2
enddo
endif
!hybr_k=hybr_k*(-2*sin(pi*k/(nz+1)))


end function hybr_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind=idp) function hybr_k_pbc(k2,k,nz)
integer, intent(in)::k,k2,nz
integer::i

hybr_k_pbc=0

!do i=1,nz
!hybr_k=hybr_k+sin(pi*i*k2/(nz+1))*cos(pi*i*k/(nz+1))
!hybr_k=hybr_k+exp(-2*pi*i*k2/(nz))*((pi*(i-1)*k/(nz+1))-sin(pi*(i+1)*k/(nz+1)))/2
!enddo

!hybr_k=hybr_k*(-2*sin(pi*k/(nz+1)))


!if (mod(k2-k-1+nz,nz)==0) hybr_k_pbc=+0.5
!if (mod(k2-k+1+nz,nz)==0) hybr_k_pbc=-0.5

if (k2==k) hybr_k_pbc=-(0,1)*sin(2*pi*k/nz)
!if (k2==k) hybr_k_pbc=-(0,1)


!print *, hybr_k_pbc

!do i=0,nz-1
!hybr_k_pbc=hybr_k_pbc+exp(-2*pi*i*k2*(0,1)/(nz)) * (exp(2*pi*(i-1)*k*(0,1)/(nz))-exp(2*pi*(i+1)*k*(0,1)/(nz))) / (2*nz)	!1/nz is the normalization
!enddo


end function hybr_k_pbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

subroutine do_dissipation
real(kind=idp),allocatable:: e1(:),e0(:)
complex(kind=idpc),allocatable::w1(:,:),w0(:,:),h_psi0(:),h_psi1(:),V_psi_tot(:,:),Ham0(:,:),Ham1(:,:),V(:,:),h_psi_tot0(:,:),h_psi_tot1(:,:),psi0vpsi0(:),psi1vpsi1(:),Vw1(:,:),Vw0(:,:)
integer,allocatable::index_psi_01(:,:),index_psi_01_tot(:,:,:)
real(kind=idp)::occup0,occup1,mu_kr0
integer::index_kr


print *, "Dissipation"


print *, "H,H0,V"

!H,H0,V
print *, nk_r

if(.not. allocated(e1)) allocate(e1(0:total_size-1))
if(.not. allocated(e0)) allocate(e0(0:total_size-1))
if(.not. allocated(w1)) allocate(w1(0:total_size-1,0:total_size-1))
if(.not. allocated(w0)) allocate(w0(0:total_size-1,0:total_size-1))


allocate(b_kr_0(0:n_sites-1))
allocate(lambda_kr_0(0:n_sites-1))
allocate(Ham1(0:total_size-1,0:total_size-1))
allocate(Ham0(0:total_size-1,0:total_size-1))
!allocate(Ham1k(0:total_size-1,0:total_size-1,0:nk_r-1))
!allocate(Ham0k(0:total_size-1,0:total_size-1,0:nk_r-1))
allocate(psi0vpsi0(0:nk_r-1))
allocate(Vw1(0:total_size-1,0:total_size-1))
allocate(Vw0(0:total_size-1,0:total_size-1))
allocate(V(0:total_size-1,0:total_size-1))
allocate(psi1vpsi1(0:nk_r-1))
allocate(index_psi_01_tot(0:total_size-1,n_max_link,3))
allocate(h_psi_tot0(0:total_size-1,n_max_link))
allocate(h_psi_tot1(0:total_size-1,n_max_link))
allocate(h_psi_tot0_k(0:total_size-1,n_max_link,0:nk_r-1))
allocate(h_psi_tot1_k(0:total_size-1,n_max_link,0:nk_r-1))
allocate(V_psi_tot (0:total_size-1,n_max_link))
allocate(index_psi_01(n_max_link,8))
allocate(h_psi0(n_max_link))
allocate(h_psi1(n_max_link))
allocate(ef_site0(0:n_sites-1,0:l_orbital_f-1))
allocate(ec_site0(0:n_sites-1,0:l_orbital-1))
allocate(ef0_site(0:n_sites-1,0:l_orbital_f-1))

psi0vpsi0=0
psi1vpsi1=0
ef0_site=0

!SCF params for H0 taken from kr calc
do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
   ind=index_xyz(ix,iy,iz,nx,ny)
   index_kr=index_xyz(mod(ix,nx_kr),mod(iy,ny_kr),iz,nx_kr,ny_kr)
   b_kr_0(ind)=b_kr(index_kr)
   lambda_kr_0(ind)=lambda_kr(index_kr)
   ef_site0(ind,:)=ef_site_kr(index_kr,:)
   ec_site0(ind,:)=ec_site_kr(index_kr,:)
   enddo
  enddo
 enddo
mu_kr0=mu_kr


if (abs(mu-mu_kr0)>1e-3) then
print *, "WARNING different mu", mu, mu_kr0
stop
endif





!loop over k_r
do ik=0,nk_r-1
!H
!this is for psi_0 and psi_1
  write (*,"(a,i5,4f10.5)") "k point", ik ,k_vecs_r(ik,1),k_vecs_r(ik,2),k_vecs_r(ik,3),k_vecs_r(ik,4)
  if (      lv_nscf) call build_ham_r(ik,Ham1,b_kr_0,lambda_kr_0,mu_kr0,ec_site,ef_site)		!nscf:uses scf from kr, not from r
  if (.not. lv_nscf) call build_ham_r(ik,Ham1,b,lambda,mu,ec_site,ef_site)
                     call build_ham_r(ik,Ham0,b_kr_0,lambda_kr_0,mu_kr0,ec_site0,ef_site0)
!  Ham0k(:,:,ik)=Ham0
!  Ham1k(:,:,ik)=Ham1
  print *, "diag_ham"
  call diagonalize_ham(Ham1,e1,w1)
  call diagonalize_ham(Ham0,e0,w0)
  print *, "OK diag_ham"


!this is for V (the time-evolution operator is not the Hamiltonian, due to the dtheta_dt term!)
  if (      lv_nscf) call build_ham_r(ik,Ham1,b_kr_0,zero,mu_kr0,ec_site,ef_site)		!nscf:uses scf from kr, not from r
  if (.not. lv_nscf) call build_ham_r(ik,Ham1,b     ,zero,mu    ,ec_site ,ef_site)
                     call build_ham_r(ik,Ham0,b_kr_0,zero,mu_kr0,ec_site0,ef_site0)
  V=Ham1-Ham0



h_psi0=0
h_psi1=0

print *, "<psi|V|psi>"

 call matmul1(total_size,V,w1,Vw1)
 call matmul1(total_size,V,w0,Vw0)


do i=0,total_size-1	!band
  occup0=1/ (exp((e0(i))/T)+1)*k_vecs_r(ik,4)/real(nk_r_tot)	!includes weight of k point mu is already in e0,e1
  occup1=1/ (exp((e1(i))/T)+1)*k_vecs_r(ik,4)/real(nk_r_tot)
  do j=0,total_size-1

    psi0Vpsi0(ik)=psi0Vpsi0(ik)+conjg(w0(j,i))*Vw0(j,i)*occup0
    psi1Vpsi1(ik)=psi1Vpsi1(ik)+conjg(w1(j,i))*Vw1(j,i)*occup1

enddo
enddo



print *, "psi0Vpsi0= ",psi0Vpsi0(ik)
print *, "psi1Vpsi1= ",psi1Vpsi1(ik)
print *, "diff= ",psi0Vpsi0(ik)-psi1Vpsi1(ik)

enddo	!ik

print *, "E_diss_tot= ",sum(psi0Vpsi0)-sum(psi1Vpsi1)

OPEN(unit=65,file=trim(label)//'/dissipation',status='unknown')

 write(65, "(3f15.6)") real(sum(psi0Vpsi0))-real(sum(psi1Vpsi1)),real(sum(psi0Vpsi0)),real(sum(psi1Vpsi1))

 close(unit=65)



end subroutine do_dissipation

!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp) function dpdr(r)
real(kind=idp),intent(in):: r

if (scf_type=='ga') dpdr=-4*r/(2-r**2)**2
if (scf_type=='sb') dpdr=-2*r

end function dpdr
!!!!!!!!
real(kind=idp) function drdp_r(r)
real(kind=idp),intent(in):: r

drdp_r=1/dpdr(r)

end function drdp_r
!!!!!!!!

real(kind=idp) function p_r(r)
real(kind=idp),intent(in):: r

if (scf_type=='ga') p_r=(1-r**2)/(1-r**2/2)
if (scf_type=='sb') p_r=(1-r**2)

end function p_r
!!!!!!!!
real(kind=idp) function r_p(p)
real(kind=idp),intent(in):: p

if (scf_type=='ga') r_p=sqrt((1-p)/(1-p/2))
if (scf_type=='sb') r_p=sqrt(1-p)

end function r_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_evolution
real ( kind = idp ), allocatable::theta_t(:),P_t(:),y(:),yp(:),b_t(:),nf_t(:),thresholds(:),dP_dt(:),dtheta_dt(:),e00(:),nc_t(:),Ec_t(:),Ef_t(:),V_t(:),V_t2(:) !,yp_f(:)
complex ( kind = idpc ), allocatable::psi_t(:,:,:),Ham0(:,:),w0(:,:)
  integer ( kind = 4 ) flag,i_step
  real ( kind = idp ):: relerr,abserr,t,t_out,t_start,E0,E1,n_tot,delta_P,n_tot_init,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot
type(rk_comm_real_1d)::comm
real(kind=idp)::tolerance
 CHARACTER(2) :: ax, ay, az
 CHARACTER(6) :: at
integer:: n_step_period


print *, "Time evolution"
print *, "Initial variables"
  neqn = 2*nx*ny*nz + 2*nk_r*total_size**2
print *, "neqn= ", neqn

 call system('mkdir ' //trim(label)//"/time_ev")

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )
  flag = 1
  t_start = 0.0d0
  t = 0.0D+00
  t_out = 0.0D+00

allocate(y(neqn))
allocate(yp(neqn))
allocate(thresholds(neqn))
!allocate(yp_f(neqn))
allocate(theta_t(0:nx*ny*nz-1))
allocate(dtheta_dt(0:nx*ny*nz-1))
allocate(P_t(0:nx*ny*nz-1))
allocate(dP_dt(0:nx*ny*nz-1))
allocate(b_t(0:nx*ny*nz-1))
allocate(psi_t(0:total_size-1,0:total_size-1,0:nk_r-1))
allocate(nff_imag(0:nx*ny*nz-1))
allocate(nfc_imag(0:nx*ny*nz-1))
allocate(nf_t(0:nx*ny*nz-1))
allocate(nc_t(0:nx*ny*nz-1))
allocate(Ef_t(0:nx*ny*nz-1))
allocate(Ec_t(0:nx*ny*nz-1))
allocate(V_t(0:nx*ny*nz-1))
allocate(V_t2(0:nx*ny*nz-1))
allocate(H_t(0:n_step))
allocate(Ham0(0:total_size-1,0:total_size-1))
allocate(w0(0:total_size-1,0:total_size-1))
allocate(e00(0:total_size-1))


tolerance=tolerance_rksuite
thresholds=threshold_rksuite



!H_0 vs H_1
delta_time=(t_stop-t_start)/n_step
n_step_period=int(period_time_ev/delta_time)
H_t=0
i=1
j=0
do while (i<n_step+1) 
 do j=0,n_step_period-1
  if (i+j<n_step+1) H_t(i+j)=1
 enddo
i=i+2*n_step_period
enddo


!do i=1,n_step
!print *, i, i*delta_time, H_t(i)
!enddo

print *, "Initial conditions"

!initial conditions
!initial wfc, eigenstate of H0 with lambda and ef
do ik=0,nk_r-1
 call build_ham_r(ik,Ham0,b_kr_0,lambda_kr_0,mu_kr,ec_site0,ef_site0)
  call diagonalize_ham(Ham0,e00,w0)
do i=0,total_size-1	!band
  occup=1/ (exp((e00(i))/T)+1)*k_vecs_r(ik,4)/real(nk_r_tot)	!includes weight of k point; mu is already in H
  psi_t(:,i,ik)=w0(:,i)*sqrt(occup)
enddo
enddo
!initial P and theta
do ind=0,nx*ny*nz-1
P_t(ind)=p_r(b_kr_0(ind))		!I have b=R=R(P) at t=0
theta_t(ind)=0.				!theta is zero at t=0
enddo



 call Pthetapsi_to_y(P_t,theta_t,psi_t,y)		!packs initial conditions into y

!if (time_integr=='rkf') 
 call r8_f2 ( t, y, yp )				!computes derivative at t=0 given y at t=0
 call y_to_Ptheta(dP_dt,dtheta_dt,yp)		!unpacks initial derivatives
if (time_integr=='rks') call setup(comm,t_start,y,t_stop,tolerance,thresholds)		!setup for rksuite

 call compute_E (n_tot_init,nf_t,nc_t,Ec_t,Ef_t,V_t,V_t2,E0,E1,P_t,b_kr_0,psi_t,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot)	!initial energy and occupation



delta_p=0
do ind=0,nx*ny*nz-1
delta_p=max(delta_p,abs(nf_t(ind)-P_t(ind)))
enddo
print *, "delta_p= ", delta_p
print *, "n_tot_init= ", n_tot_init


OPEN(unit=562,file=trim(label)//'/time_ev/E0E1',status='unknown')
write ( 562, '(i4,1000f14.6)' ) 0, t,n_tot_init, E0,E1,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot,P_t,b_kr_0,theta_t


!!!!t xyz files
WRITE(UNIT=at,FMT='(f0.2)') t
OPEN(unit=1000,file=trim(label)//'/time_ev/t'//at,status='unknown')
OPEN(unit=1001,file=trim(label)//'/time_ev/b_rt'//at,status='unknown')
write(1000,"(3a4,100a15)") "# ix","iy","iz", "P_t", "b_kr_0","theta_t","dP_dt", "dtheta_dt", "ef-dtheta_dt","nf_t","nc_t","Ec_t","Ef_t","V_t","V_t2"
do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
 WRITE(UNIT=ax,FMT='(I2.2)') ix
 WRITE(UNIT=ay,FMT='(I2.2)') iy
 WRITE(UNIT=az,FMT='(I2.2)') iz
 ind=index_xyz(ix,iy,iz,nx,ny)
 write(1000,"(3i4,100f15.10)") ix,iy,iz, P_t(ind), b_kr_0(ind),theta_t(ind),dP_dt(ind), dtheta_dt(ind),ef_site(ind,0)-dtheta_dt(ind),nf_t(ind),nc_t(ind),Ec_t(ind),Ef_t(ind),V_t(ind),V_t2(ind)
 write(1001,"(3i4,100f15.10)") ix,iy,iz, b_kr_0(ind), ef_site(ind,0)-dtheta_dt(ind)
 OPEN(unit=2000+ind,file=trim(label)//'/time_ev/nx'//ax//'_ny'//ay//'_nz'//az,status='unknown')
 write ( 2000+ind, '(100a14)' ) "#            t","P_t","b_t","theta_t","dP_dt", "dtheta_dt", "ef-dtheta_dt","nf_t","nc_t","Ec_t","Ef_t","V_t","V_t2"
 write ( 2000+ind, '(1000f14.6)' ) t,P_t(ind),b_kr_0(ind),theta_t(ind),dP_dt(ind), dtheta_dt(ind),ef_site(ind,0)-dtheta_dt(ind),nf_t(ind),nc_t(ind),Ec_t(ind),Ef_t(ind),V_t(ind),V_t2(ind)
enddo
enddo
enddo
 close(unit=1000)
 close(unit=1001)


print *, "Starting t loop"

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = idp ) * t_start &
        + real (          i_step - 1, kind = idp ) * t_stop ) &
        / real ( n_step,              kind = idp )

    t_out = ( real ( n_step - i_step, kind = idp ) * t_start &
            + real (          i_step, kind = idp ) * t_stop ) &
            / real ( n_step,          kind = idp )

if (time_integr=='rkf')    call r8_rkf45 ( r8_f2, neqn, y, yp, t, t_out, relerr, abserr, flag )
if (time_integr=='rks')    call range_integrate(comm,yp_f,t_out,t,y,yp,flag)

    call y_to_Pthetapsi(P_t,theta_t,psi_t,y)
    do ind=0,nx*ny*nz-1
     b_t(ind)=r_p(P_t(ind))		!R=b from P
    enddo

    call compute_E (n_tot,nf_t,nc_t,Ec_t,Ef_t,V_t,V_t2,E0,E1,P_t,b_t,psi_t,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot)

!check that P=nf (ga) or b**2+nf (sb) is conserved by the time evolution
    delta_p=0
    do ind=0,nx*ny*nz-1
    if (scf_type=='ga')    delta_p=max(delta_p,abs(nf_t(ind)-P_t(ind)))
    if (scf_type=='sb')    delta_p=max(delta_p,abs(nf_t(ind)+b_t(ind)**2-1))
    enddo
    if (scf_type=='ga' .and. delta_p>1e-4) print *, "Warning |P-nf|=  ", delta_p
    if (scf_type=='sb' .and. delta_p>1e-4) print *, "Warning |b**2+nf-1=|  ", delta_p
!check that total number of particles is conserved by the time evolution
    if (abs(n_tot-n_tot_init)>1e-3) print *, "Warning n_tot=  ", n_tot


    write ( *, '(2i4,2x,4f14.6)' ) flag, H_t(int(t/delta_time)), t,E0,E1, b_t(((nx-1)*(nx+1))/2)
    write ( 562, '(i4,1000f14.6)' ) H_t(int(t/delta_time)), t, n_tot,E0,E1,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot,P_t,b_t,theta_t

!!!!t xyz files
 call y_to_Ptheta(dP_dt,dtheta_dt,yp)		!unpacks derivatives
WRITE(UNIT=at,FMT='(f0.2)') t
OPEN(unit=1000,file=trim(label)//'/time_ev/t'//at,status='unknown')
OPEN(unit=1001,file=trim(label)//'/time_ev/b_rt'//at,status='unknown')
write(1000,"(3a4,100a15)") "# ix","iy","iz", "P_t", "b_t","theta_t","dP_dt", "dtheta_dt", "ef-dtheta_dt","nf_t","nc_t","Ec_t","Ef_t","V_t","V_t2"
do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
 ind=index_xyz(ix,iy,iz,nx,ny)
 write(1000,"(3i4,100f15.10)") ix,iy,iz, P_t(ind), b_t(ind),theta_t(ind),dP_dt(ind),dtheta_dt(ind),ef_site(ind,0)- dtheta_dt(ind),nf_t(ind),nc_t(ind),Ec_t(ind),Ef_t(ind),V_t(ind),V_t2(ind)
 write(1001,"(3i4,100f15.10)") ix,iy,iz, b_t(ind), ef_site(ind,0)- dtheta_dt(ind)
 write ( 2000+ind, '(1000f14.6)' ) t,P_t(ind),b_t(ind),theta_t(ind),dP_dt(ind),dtheta_dt(ind),ef_site(ind,0)- dtheta_dt(ind),nf_t(ind),nc_t(ind),Ec_t(ind),Ef_t(ind),V_t(ind),V_t2(ind)
enddo
enddo
enddo
 close(unit=1000)
 close(unit=1001)


  end do

 close(unit=562)
do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
 ind=index_xyz(ix,iy,iz,nx,ny)
 close(unit=2000+ind)
enddo
enddo
enddo


deallocate(y)
deallocate(yp)
deallocate(thresholds)
deallocate(theta_t)
deallocate(P_t)
deallocate(dtheta_dt)
deallocate(dP_dt)
deallocate(b_t)
deallocate(psi_t)
deallocate(nff_imag)
deallocate(nfc_imag)
deallocate(nf_t)
deallocate(nc_t)
deallocate(Ef_t)
deallocate(Ec_t)
deallocate(V_t)
deallocate(V_t2)
!deallocate(yp_f)

end subroutine time_evolution
!!!!!!!
subroutine r8_f2 ( t, y, yp )
!! R8_F2 evaluates the derivative for the ODE.
real ( kind = idp ):: t,y(neqn),yp(neqn)
complex( kind = idpc ), allocatable::Ham(:,:),Hampsi(:,:),psi_t(:,:,:),Hampsi_tot(:,:,:)
real ( kind = idp ), allocatable:: P_t(:),theta_t(:),dPdt(:),dthetadt(:),b_t(:)

time_r8f2_0=secnds(0.0)
allocate(P_t(0:n_sites-1))
allocate(b_t(0:n_sites-1))
allocate(theta_t(0:n_sites-1))
allocate(dPdt(0:n_sites-1))
allocate(dthetadt(0:n_sites-1))
allocate(psi_t(0:total_size-1,0:total_size-1,0:nk_r-1))

 call y_to_Pthetapsi(P_t,theta_t,psi_t,y)	!unpacks P theta and psi from y

do ind=0,nx*ny*nz-1		!finds b=R from P
 b_t(ind)=R_p(P_t(ind))
enddo

 call compute_nf_r_imag (nf,nc,nfc,nff,nfc_imag, nff_imag,b_t,psi_t)	!computes <f+c> etc from psi_t and b

do ind=0,nx*ny*nz-1
dPdt(ind)=-(nfc_imag(ind)+nff_imag(ind))*b_t(ind)
if (H_t(int(t/delta_time))==1) dthetadt(ind)=(nfc(ind)+nff(ind))*drdP_r(b_t(ind))+ef_site(ind,0)		!ef_site has as second index the number of f orbital, but here I only consider 1 f orbital
if (H_t(int(t/delta_time))==0) dthetadt(ind)=(nfc(ind)+nff(ind))*drdP_r(b_t(ind))+ef_site0(ind,0)		!ef_site has as second index the number of f orbital, but here I only consider 1 f orbital
enddo

allocate(Ham(0:total_size-1,0:total_size-1))
allocate(Hampsi(0:total_size-1,0:total_size-1))
allocate(Hampsi_tot(0:total_size-1,0:total_size-1,0:nk_r-1))

!dpsi_dt
do ik=0,nk_r-1
 call build_ham_r(ik,Ham,b_t,-dthetadt,mu,ec_site,ef0_site)		!I do not have lambda in H(t) or ef nf, but an additional -dtheta/dt n   !ef0_site=0
 call matmul1(total_size,Ham,psi_t(:,:,ik),Hampsi)
Hampsi_tot(:,:,ik)=(0,1.)*Hampsi
enddo


 call Pthetapsi_to_y(dPdt,dthetadt,Hampsi_tot,yp)		!packs derivatives into yp

deallocate(Ham)
deallocate(Hampsi)
deallocate(Hampsi_tot)
deallocate(P_t)
deallocate(b_t)
deallocate(theta_t)
deallocate(dPdt)
deallocate(dthetadt)
deallocate(psi_t)

!print *, "ok r8_f2"
time_r8f2=time_r8f2+secnds(0.0)-time_r8f2_0

end subroutine r8_f2

!!!!!!
!wrapper for r8_f2: function that returns yp
function yp_f2( t, y) !result(res)
!! R8_F2 evaluates the derivative for the ODE.
  real ( kind = idp ) t
  real ( kind = idp ),dimension(:),intent(in):: y
  real ( kind = idp ),allocatable:: res(:)
  real ( kind = idp ),dimension(size(y)):: yp_f2

allocate(res(size(y)))
 call r8_f2(t,y,res)
yp_f2=res
deallocate(res)

end function yp_f2

!!!!!!!!!!!!!!!!!!!
function yp_f( t, y)
!! R8_F2 evaluates the derivative for the ODE.
  real ( kind = idp ) t
  real ( kind = idp ),dimension(:),intent(in):: y
  real ( kind = idp )::yp_f(size(y))

 call r8_f2(t,y,yp_f)


end function yp_f

!!!!!!!!!!!!!!!!!!!
subroutine compute_nf_r_imag (nf,nc,nfc,nff,nfc_imag, nff_imag,b,w_r_tot)
complex(KIND=idpc) ,intent(in)	:: w_r_tot(0:total_size-1,0:total_size-1, 0:nk_r-1) !, h_psi_tot(0:total_size-1,n_max_link,0:nk_r-1)
real(KIND=idp) ,intent(in)	:: b(0:n_sites-1)
real(kind=idp) ,intent(out)	:: nf(0:n_sites-1),nc(0:n_sites-1),nfc(0:n_sites-1),nff(0:n_sites-1),nfc_imag(0:n_sites-1),nff_imag(0:n_sites-1)
real(kind=idp)			:: occup
integer				::i,j,ix,iy,iz,icf,is,o1
integer				::ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8
real(kind=idp) ,allocatable	:: nf_corr(:),nc_corr(:),nfc_corr(:),nff_corr(:),nfc_imag_corr(:),nff_imag_corr(:)

time_nfimag_0=secnds(0.0)

!h_psi_tot is the non-diagonal part of the Hamiltonian, that I assume equal for H0 and H1
 
nf=0
nc=0
nfc=0	!V included
nff=0	!b*tf included
nfc_imag=0	!V included
nff_imag=0	!b*tf included

 do ik=0, nk_r-1  ! k points
!  print *, "k=", ik
   do i=0, total_size-1	!eigenvalues
!    occup=1/ (exp((e_r_tot(i,ik)-mu)/T)+1)*k_vecs_r(ik,4)/real(nk_r_tot) !I already put this in w_r_tot
     occup=1
    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z
       do icf=0,1	!c/f
        do is=0,1	!up/down 
         do o1=0,l_orb(icf)-1

	ind=index_xyz(ix,iy,iz,nx,ny)					!site
	indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz)
	
!occupation
if (icf==1) nf (ind)= nf (ind) + (abs(w_r_tot(indp,i,ik)))**2 * occup
if (icf==0) nc (ind)= nc (ind) + (abs(w_r_tot(indp,i,ik)))**2 * occup	

!!!nfc,nff

     do j=1,n_max_link
      ind2 =index_psi_tot(indp,j,1)
      ind2p=index_psi_tot(indp,j,2)
      icf2 =index_psi_tot(indp,j,3)
    
  
if (icf==1 .and. icf2==0) nfc(ind)=nfc(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   *occup
if (icf==1 .and. icf2==1) nff(ind)=nff(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		   	 	            *b(ind2) &
		  	  		   *occup
if (icf==1 .and. icf2==0) nfc_imag(ind)=nfc_imag(ind)+aimag( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 -conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   *occup
if (icf==1 .and. icf2==1) nff_imag(ind)=nff_imag(ind)+aimag( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 -conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		   	 	            *b(ind2) &
		  	  		   *occup


	   enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo !ik
   
!correct occupations for reduced set of k points

if (lkred_r .and. nx==ny .and. nx>1 .and. lysameasx .and. .not. ldisorder_r) then

allocate(nc_corr(0:n_sites-1))
allocate(nf_corr(0:n_sites-1))
allocate(nfc_corr(0:n_sites-1))
allocate(nff_corr(0:n_sites-1))
allocate(nfc_imag_corr(0:n_sites-1))
allocate(nff_imag_corr(0:n_sites-1))

    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z

	ind=ix + iy*nx + iz*nx*ny		!ij
	ind_2=(nx-1-ix) + iy*nx + iz*nx*ny	!-ij
	ind_3=ix + (ny-1-iy)*nx + iz*nx*ny		!i-j
	ind_4=(nx-1-ix) + (ny-1-iy)*nx + iz*nx*ny		!-i-j
	ind_5=iy + ix*nx + iz*nx*ny		!ji
	ind_6=(ny-1-iy) + ix*nx + iz*nx*ny		!-ji
	ind_7=iy + (nx-1-ix)*nx + iz*nx*ny		!j-i
	ind_8=(ny-1-iy) + (nx-1-ix)*nx + iz*nx*ny		!-j-i


nf_corr(ind)=(nf(ind)+nf(ind_2)+nf(ind_3)+nf(ind_4)+nf(ind_5)+nf(ind_6)+nf(ind_7)+nf(ind_8))/8
nc_corr(ind)=(nc(ind)+nc(ind_2)+nc(ind_3)+nc(ind_4)+nc(ind_5)+nc(ind_6)+nc(ind_7)+nc(ind_8))/8
nfc_corr(ind)=(nfc(ind)+nfc(ind_2)+nfc(ind_3)+nfc(ind_4)+nfc(ind_5)+nfc(ind_6)+nfc(ind_7)+nfc(ind_8))/8
nff_corr(ind)=(nff(ind)+nff(ind_2)+nff(ind_3)+nff(ind_4)+nff(ind_5)+nff(ind_6)+nff(ind_7)+nff(ind_8))/8
nfc_imag_corr(ind)=(nfc_imag(ind)+nfc_imag(ind_2)+nfc_imag(ind_3)+nfc_imag(ind_4)+nfc_imag(ind_5)+nfc_imag(ind_6)+nfc_imag(ind_7)+nfc_imag(ind_8))/8
nff_imag_corr(ind)=(nff_imag(ind)+nff_imag(ind_2)+nff_imag(ind_3)+nff_imag(ind_4)+nff_imag(ind_5)+nff_imag(ind_6)+nff_imag(ind_7)+nff_imag(ind_8))/8

enddo
enddo
enddo

nf=nf_corr
nc=nc_corr
nfc=nfc_corr
nff=nff_corr
nfc_imag=nfc_imag_corr
nff_imag=nff_imag_corr


deallocate(nc_corr)
deallocate(nf_corr)
deallocate(nfc_corr)
deallocate(nff_corr)
deallocate(nfc_imag_corr)
deallocate(nff_imag_corr)
endif

time_nfimag=time_nfimag+secnds(0.0)-time_nfimag_0

end subroutine compute_nf_r_imag
!!!!!!!!!!!!!!!!!!!
subroutine Pthetapsi_to_y (P,theta,psi,y)
real(kind=idp), intent(in):: P(0:nx*ny*nz-1), theta(0:nx*ny*nz-1)
complex(kind=idpc), intent(in):: psi(0:total_size-1,0:total_size-1,0:nk_r-1)
real(kind=idp), intent(out):: y(neqn)
logical,allocatable:: check(:)
integer::ind_p, ind_t,ind_psir,ind_psii

allocate(check(neqn))
 check=.false.

do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
   ind=index_xyz(ix,iy,iz,nx,ny)
   ind_p=ind+1
   y(ind_p)=P(ind)
   check(ind_p)=.true.
   ind_t=ind_p+nx*ny*nz
   y(ind_t)=theta(ind)
   check(ind_t)=.true.
   enddo
  enddo
 enddo

!psi
do i=0,total_size-1	
  do n=0,total_size-1	!eigenvalue
    do ik=0,nk_r-1
ind_psir=2*nx*ny*nz+i+total_size*n+total_size**2*ik+1
ind_psii=ind_psir+total_size**2*nk_r
y(ind_psir)=real(psi(i,n,ik))
y(ind_psii)=aimag(psi(i,n,ik))
 check(ind_psir)=.true.
 check(ind_psii)=.true.
enddo
enddo
enddo

do i=1,neqn
if(.not. (check(i))) print *, "pthetapsi_to_y",  i, check(i)
enddo

deallocate(check)


end subroutine pthetapsi_to_y

!!!!!!!!!!!!!!!
subroutine y_to_Pthetapsi (P,theta,psi,y)
real(kind=idp), intent(out):: P(0:nx*ny*nz-1), theta(0:nx*ny*nz-1)
complex(kind=idp), intent(out):: psi(0:total_size-1,0:total_size-1,0:nk_r-1)
real(kind=idp), intent(in):: y(neqn)
logical,allocatable:: check(:)
integer::ind_p, ind_t,ind_psir,ind_psii,total_ind

allocate(check(neqn))
 check=.false.
total_ind=0

do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
   ind=index_xyz(ix,iy,iz,nx,ny)
   ind_p=ind+1
   ind_t=ind_p+nx*ny*nz
   P(ind)=y(ind_p)
   theta(ind)=y(ind_t)
   check(ind_p)=.true.
   check(ind_t)=.true.
total_ind=total_ind+1
   enddo
  enddo
 enddo

!psi
do i=0,total_size-1
  do n=0,total_size-1
    do ik=0,nk_r-1

ind_psir=2*nx*ny*nz+i+total_size*n+total_size**2*ik+1
ind_psii=ind_psir+total_size**2*nk_r
psi(i,n,ik)=y(ind_psir)+(0,1.)*y(ind_psii)
 check(ind_psir)=.true.
 check(ind_psii)=.true.
total_ind=total_ind+1
enddo
enddo
enddo


do i=1,neqn
if(.not. (check(i))) print *, "y_to_pthetapsi",  i, check(i)
enddo
if(2*total_ind .ne. neqn) print *, total_ind, neqn

deallocate(check)

end subroutine y_to_Pthetapsi

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine y_to_Ptheta (P,theta,y)
real(kind=idp), intent(out):: P(0:nx*ny*nz-1), theta(0:nx*ny*nz-1)
real(kind=idp), intent(in):: y(neqn)
integer::ind,ind_t,ind_p

do ix=0,nx-1
  do iy=0,ny-1
   do iz=0,nz-1
   ind=index_xyz(ix,iy,iz,nx,ny)
   ind_p=ind+1
   ind_t=ind_p+nx*ny*nz
   P(ind)=y(ind_p)
   theta(ind)=y(ind_t)
   enddo
  enddo
 enddo

end subroutine y_to_Ptheta

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_E (n_tot,nf_t,nc_t,Ec_t,Ef_t,V_t,V_t2,E0,E1,P_t,b_t,w_r_tot,Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot)
complex(KIND=idpc) ,intent(in)	:: w_r_tot(0:total_size-1,0:total_size-1, 0:nk_r-1) !, h_psi_tot(0:total_size-1,n_max_link,0:nk_r-1)
real(KIND=idp) ,intent(in)	:: P_t(0:n_sites-1),b_t(0:n_sites-1)
real(kind=idp) ,intent(out)	:: E0,E1,n_tot,nf_t(0:nx*ny*nz-1),nc_t(0:nx*ny*nz-1),Ef_t(0:nx*ny*nz-1),V_t(0:nx*ny*nz-1),Ec_t(0:nx*ny*nz-1),V_t2(0:nx*ny*nz-1),Ec0_tot,Ec1_tot,Tc_tot,Ef0_tot,Ef1_tot,Tf_tot,V_tot
integer				::i,j,ix,iy,iz,icf,is,o1,ik
integer				::ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8
complex(KIND=idpc),allocatable:: w_r(:,:),Hw0(:,:),Hampsi(:,:),Hw1(:,:),H0k(:,:),H1k(:,:),Hdk(:,:)
complex(KIND=idpc),allocatable:: Ec_corr(:),Ef_corr(:),V_corr(:),V2_corr(:),nf_corr(:),nc_corr(:)

time_e_0=secnds(0.0)

allocate(Hw0(0:total_size-1,0:total_size-1))
allocate(Hw1(0:total_size-1,0:total_size-1))
allocate(H0k(0:total_size-1,0:total_size-1))
allocate(H1k(0:total_size-1,0:total_size-1))
allocate(Hdk(0:total_size-1,0:total_size-1))
allocate(w_r(0:total_size-1,0:total_size-1))

E0=0
E1=0
n_tot=0
nf_t=0
nc_t=0
Ec_t=0
Ef_t=0
V_t=0
V_t2=0

Tc_tot=0
Ec0_tot=0
Ec1_tot=0
Tf_tot=0
Ef0_tot=0
Ef1_tot=0
V_tot=0


do ik=0,nk_r-1
w_r=w_r_tot(:,:,ik)		!this already includes occupations

 call build_ham_r(ik,H0k,b_t,zero,mu_kr,ec_site0,ef_site0)		!ef0_site is zero ef_site0 is the kr value (H0)
 call build_ham_r(ik,H1k,b_t,zero,mu_kr,ec_site ,ef_site)		!ef_site is the H1 value
! call build_ham_r(ik,Hdk,b_t,zero,mu_kr,ec_site ,ef0_site)		!ef_site is the H1 value
 call matmul1(total_size,H0k,w_r,Hw0)
 call matmul1(total_size,H1k,w_r,Hw1)

!total energy
 do i=0,total_size-1	!band
 do j=0,total_size-1
  E0=E0+conjg(w_r(j,i))*Hw0(j,i)
  E1=E1+conjg(w_r(j,i))*Hw1(j,i)
  n_tot=n_tot+conjg(w_r(j,i))*w_r(j,i) 
 enddo
 enddo



   do i=0, total_size-1	!eigenvalues
    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z
       do icf=0,1	!c/f
        do is=0,1	!up/down 
         do o1=0,l_orb(icf)-1

	ind=index_xyz(ix,iy,iz,nx,ny)					!site
	indp=index_xyzcso(ix,iy,iz,icf,is,o1,nx, ny, nz)

if (icf==0) nc_t (ind)= nc_t (ind) + (abs(w_r_tot(indp,i,ik)))**2
if (icf==1) nf_t (ind)= nf_t (ind) + (abs(w_r_tot(indp,i,ik)))**2

!global contributions to energy
if (icf==0) Ec0_tot = Ec0_tot+ (abs(w_r_tot(indp,i,ik)))**2*(ec_site0(ind,0)-mu_kr)
if (icf==0) Ec1_tot = Ec1_tot+ (abs(w_r_tot(indp,i,ik)))**2*(ec_site(ind,0)-mu_kr)
if (icf==1) Ef0_tot = Ef0_tot+ (abs(w_r_tot(indp,i,ik)))**2*(ef_site0(ind,0)-mu_kr)
if (icf==1) Ef1_tot = Ef1_tot+ (abs(w_r_tot(indp,i,ik)))**2*(ef_site(ind,0)-mu_kr)


     do j=1,n_max_link
      ind2 =index_psi_tot(indp,j,1)
      ind2p=index_psi_tot(indp,j,2)
      icf2 =index_psi_tot(indp,j,3)
    
!local contributions to energy
if ((icf==1 .and. icf2==0)) &		
			  V_t(ind)=V_t(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   *b_t(ind)
if ((icf==0 .and. icf2==1)) &		
			  V_t2(ind)=V_t2(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   *b_t(ind2)

if (icf==1 .and. icf2==1) Ef_t(ind)=Ef_t(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		   	 	            *b_t(ind2)*b_t(ind) &
		  	  		   /2

if (icf==0 .and. icf2==0) Ec_t(ind)=Ec_t(ind)+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  		   /2



!global contributions to energy
!/2 to avoid double counting
if ((icf==1 .and. icf2==0)) 	  V_tot=V_tot+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	  +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		  	  	  *b_t(ind)

if (icf==0 .and. icf2==0) 	Tc_tot=Tc_tot+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik))/2


if (icf==1 .and. icf2==1) Tf_tot=Tf_tot+real( conjg(w_r_tot(indp,i,ik)) *h_psi_tot(indp,j,ik)*w_r_tot(ind2p,i,ik)&
		   	 	                 +conjg(w_r_tot(ind2p,i,ik))*conjg(h_psi_tot(indp,j,ik))*w_r_tot(indp,i,ik)) &
		   	 	            *b_t(ind2)*b_t(ind)/2



enddo	!Links

enddo
enddo
enddo
enddo
enddo
enddo
enddo




enddo	!ik

!P_i*ef_fi
!do ind=0,nx*ny*nz-1
!Ef_t(ind)=Ef_t(ind)+P_t(ind)*(ef_site0(ind,0)-mu)
!Ec_t(ind)=Ec_t(ind)+nc_t(ind)*(ec_site0(ind,0)-mu)
!E0=E0+P_t(ind)*ef_site0(ind,0)
!E1=E1+P_t(ind)*ef_site(ind,0)
!enddo



allocate(Ef_corr(0:n_sites-1))
allocate(Ec_corr(0:n_sites-1))
allocate(V_corr(0:n_sites-1))
allocate(V2_corr(0:n_sites-1))
allocate(nf_corr(0:n_sites-1))
allocate(nc_corr(0:n_sites-1))

    do ix=0,nx-1	!x
     do iy=0,ny-1	!y
      do iz=0,nz-1	!z

	ind=ix + iy*nx + iz*nx*ny		!ij
	ind_2=(nx-1-ix) + iy*nx + iz*nx*ny	!-ij
	ind_3=ix + (ny-1-iy)*nx + iz*nx*ny		!i-j
	ind_4=(nx-1-ix) + (ny-1-iy)*nx + iz*nx*ny		!-i-j
	ind_5=iy + ix*nx + iz*nx*ny		!ji
	ind_6=(ny-1-iy) + ix*nx + iz*nx*ny		!-ji
	ind_7=iy + (nx-1-ix)*nx + iz*nx*ny		!j-i
	ind_8=(ny-1-iy) + (nx-1-ix)*nx + iz*nx*ny		!-j-i


nc_corr(ind)=(nc_t(ind)+nc_t(ind_2)+nc_t(ind_3)+nc_t(ind_4)+nc_t(ind_5)+nc_t(ind_6)+nc_t(ind_7)+nc_t(ind_8))/8
nf_corr(ind)=(nf_t(ind)+nf_t(ind_2)+nf_t(ind_3)+nf_t(ind_4)+nf_t(ind_5)+nf_t(ind_6)+nf_t(ind_7)+nf_t(ind_8))/8
Ec_corr(ind)=(Ec_t(ind)+Ec_t(ind_2)+Ec_t(ind_3)+Ec_t(ind_4)+Ec_t(ind_5)+Ec_t(ind_6)+Ec_t(ind_7)+Ec_t(ind_8))/8
Ef_corr(ind)=(Ef_t(ind)+Ef_t(ind_2)+Ef_t(ind_3)+Ef_t(ind_4)+Ef_t(ind_5)+Ef_t(ind_6)+Ef_t(ind_7)+Ef_t(ind_8))/8
V_corr(ind)=(V_t(ind)+V_t(ind_2)+V_t(ind_3)+V_t(ind_4)+V_t(ind_5)+V_t(ind_6)+V_t(ind_7)+V_t(ind_8))/8
V2_corr(ind)=(V_t2(ind)+V_t2(ind_2)+V_t2(ind_3)+V_t2(ind_4)+V_t2(ind_5)+V_t2(ind_6)+V_t2(ind_7)+V_t2(ind_8))/8
enddo
enddo
enddo

nf_t=nf_corr
nc_t=nc_corr
Ef_t=Ef_corr
Ec_t=Ec_corr
V_t=V_corr
V_t2=V2_corr



deallocate(nc_corr)
deallocate(nf_corr)
deallocate(Ec_corr)
deallocate(Ef_corr)
deallocate(V_corr)
deallocate(V2_corr)	




deallocate(Hw0)
deallocate(Hw1)
deallocate(H0k)
deallocate(H1k)
deallocate(w_r)

time_e=time_e+secnds(0.0)-time_e_0

end subroutine compute_E



end program tki


