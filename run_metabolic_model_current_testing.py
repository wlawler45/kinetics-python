import kinetics
import matplotlib.pyplot as plt
import reaction_equations
if __name__ == "__main__":
	model = kinetics.Model(logging=False)
	# vPts1 = reaction_equations.vPts1_Reaction()
	# vPts1.parameters={'EIIAtotal': 0.0769000000,'kmPts1': 27700.0000000000,'kPts1': 66000.0000000000}
	# model.append(vPts1)
	# vPts4 = reaction_equations.vPts4_Reaction()
	# vPts4.parameters={'KPts_GLC': 0.0049000000,'KPts_EIIA': 0.0021100000,'vPts4_max': 3940.0000000000}
	# model.append(vPts4)
	# vPts4_medium = reaction_equations.vPts4_medium_Reaction()
	# vPts4_medium.parameters={'KPts_GLC': 0.0049000000,'KPts_EIIA': 0.0021100000,'vPts4_max': 3940.0000000000,'rho': 564.0000000000}
	# model.append(vPts4_medium)
	# vNonpts = reaction_equations.vNonpts_Reaction()
	# vNonpts.parameters={'KNonpts_S': 1.5500000000,'EIIAtotal': 0.0769000000,'vNonpts_max': 4000.0000000000,'KNonpts_I': 0.0100000000}
	# model.append(vNonpts)
	# vNonpts_medium = reaction_equations.vNonpts_medium_Reaction()
	# vNonpts_medium.parameters={'KNonpts_S': 1.5500000000,'rho': 564.0000000000,'EIIAtotal': 0.0769000000,'vNonpts_max': 4000.0000000000,'KNonpts_I': 0.0100000000}
	# model.append(vNonpts_medium)
	vE_Glk = reaction_equations.vE_Glk_Reaction()
	vE_Glk.parameters = {'KGlk_GLC_m': 220.0000000000, 'KGlk_ATP_m': 800.0000000000, 'KGlk_G6P_i': 15000.0000000000,
						 'kGlk_cat': 6.8333333333}
	model.append(vE_Glk)
	vE_Pgi = reaction_equations.vE_Pgi_Reaction()
	vE_Pgi.parameters = {'vPgi_max': 39000000.0000000000, 'KPgi_G6P_6pginh': 181.0000000000, 'KPgi_G6P': 2460.0000000000,
						 'KPgi_F6P_6pginh': 199.0000000000, 'KPgi_F6P': 342.0000000000, 'KPgi_eq': 1440.0000000000}
	model.append(vE_Pgi)
	vE_Pfk = reaction_equations.vE_Pfk_Reaction()
	vE_Pfk.parameters = {'KPfk_AMP_a': 10100.0000000000, 'KPfk_AMP_b': 25.5000000000, 'LPfk': 1770000.0000000000,
						 'KPfk_PEP': 1740.0000000000, 'KPfk_ADP_b': 256.0000000000, 'KPfk_ADP_a': 277000.0000000000,
						 'KPfk_ADP_c': 448.0000000000, 'KPfk_ATP_s': 161.0000000000, 'KPfk_F6P_s': 21.1000000000,
						 'kPfk_cat': 77833333.3333333284, 'nPfk': 4.0000000000}
	model.append(vE_Pfk)
	vE_Fbp = reaction_equations.vE_Fbp_Reaction()
	vE_Fbp.parameters = {'LFbp': 4410000.0000000000, 'KFbp_PEP': 488.0000000000, 'KFbp_FBP': 8.9200000000,
						 'nFbp': 4.0000000000, 'kFbp_cat': 1310000.0000000000}
	model.append(vE_Fbp)
	vE_Fba = reaction_equations.vE_Fba_Reaction()
	vE_Fba.parameters = {'KFba_eq': 372.0000000000, 'KFba_DHAP': 88.0000000000, 'kFba_cat': 11583333.3333333340,
						 'VFba_blf': 1.5400000000, 'KFba_FBP': 83.9000000000, 'KFba_GAP_inh': 600.0000000000,
						 'KFba_GAP': 154.0000000000}
	model.append(vE_Fba)
	vE_Eno = reaction_equations.vE_Eno_Reaction()
	vE_Eno.parameters = {'KEno_eq': 2999.0000000000, 'vEno_max': 0.1953150000, 'KEno_Pep': 100.0000000000,
						 'KEno_Pga2': 100.0000000000}
	model.append(vE_Eno)
	vE_Tpi = reaction_equations.vE_Tpi_Reaction()
	vE_Tpi.parameters = {'K_TPI_KmGAP': 1893.0100000000, 'K_TPI_Keq': 270.2030000000, 'K_TPI_KmDAP': 10.0000000000,
						 'v_TPI_Vmax': 0.4030716667}
	model.append(vE_Tpi)
	vE_Pgk = reaction_equations.vE_Pgk_Reaction()
	vE_Pgk.parameters = {'K_PGK_KmPGA3': 2457.2200000000, 'K_PGK_KmADPMg': 85.4160000000, 'v_PGK_Vmax': 0.2684816667,
						 'K_PGK_Keq': 99992.5000000000, 'K_PGK_KmATPMg': 3477.3700000000, 'K_PGK_KmBPG': 11.3296000000}
	model.append(vE_Pgk)
	vE_Gpm = reaction_equations.vE_Gpm_Reaction()
	vE_Gpm.parameters = {'K_GPM_KmPGA3': 115.0000000000, 'v_GPM_Vmax': 0.1832233333, 'K_GPM_Keq': 565.8180000000,
						 'K_GPM_KmPGA2': 1915.3000000000}
	model.append(vE_Gpm)
	vE_Gdh = reaction_equations.vE_Gdh_Reaction()
	vE_Gdh.parameters = {'K_GDH_Keq': 20000.0000000000, 'K_GDH_KmNADH': 3697.9700000000, 'K_GDH_KmBPG': 200.0000000000,
						 'K_GDH_KmP': 17.0000000000, 'v_GDH_Vmax': 0.1444288333, 'K_GDH_KmGAP': 2472.6500000000,
						 'K_GDH_KmNAD': 11.0454000000}
	model.append(vE_Gdh)
	vE_Pyk = reaction_equations.vE_Pyk_Reaction()
	vE_Pyk.parameters = {'LPyk': 997.0000000000, 'KPyk_ATP': 20200.0000000000, 'KPyk_ADP': 210.0000000000,
						 'nPyk': 4.0000000000, 'kPyk_cat': 13600.0000000000, 'KPyk_AMP': 255.0000000000,
						 'KPyk_FBP': 253.0000000000, 'KPyk_PEP': 310.0000000000}
	model.append(vE_Pyk)
	# vE_Pps = reaction_equations.vE_Pps_Reaction()
	# vE_Pps.parameters={'nPps': 2.0000000000,'KPps_PEP': 0.0002160000,'KPps_PYR': 0.0007130000,'kPps_cat': 249000.0000000000,'LPps': 0.0000000000}
	# model.append(vE_Pps)
	#vE_Pdh = reaction_equations.vE_Pdh_Reaction()
	#vE_Pdh.parameters={'KPdh_CoA_m': 0.0048200000,'kPdh_cat': 48300000.0000000000,'KPdh_AcCoA_m': 0.0080000000,'KPdh_i': 68.3000000000,'KPdh_NADH_m': 0.0471000000,'KPdh_NAD_m': 0.4010000000,'KPdh_PYR_m': 1.0000000000}
	#model.append(vE_Pdh)
	# vE_Pta = reaction_equations.vE_Pta_Reaction()
	# vE_Pta.parameters={'KPta_AcP_m': 0.2290000000,'KPta_CoA_i': 0.0810000000,'KPta_Pi_m': 0.6870000000,'KPta_AcP_i': 0.3200000000,'vPta_max': 5360.0000000000,'KPta_AcCoA_i': 0.2010000000,'KPta_Pi_i': 2.1100000000,'KPta_eq': 0.0284000000}
	# model.append(vE_Pta)
	# vE_Ack = reaction_equations.vE_Ack_Reaction()
	# vE_Ack.parameters={'KAck_AcP_m': 0.0463000000,'KAck_eq': 234.0000000000,'KAck_ADP_m': 0.1770000000,'KAck_ATP_m': 0.0886000000,'KAck_ACE_m': 6.1000000000,'vAck_max': 195000.0000000000}
	# model.append(vE_Ack)
	# vE_Ack_medium = reaction_equations.vE_Ack_medium_Reaction()
	# vE_Ack_medium.parameters={'KAck_AcP_m': 0.0463000000,'rho': 564.0000000000,'KAck_eq': 234.0000000000,'KAck_ADP_m': 0.1770000000,'KAck_ATP_m': 0.0886000000,'KAck_ACE_m': 6.1000000000,'vAck_max': 195000.0000000000}
	# model.append(vE_Ack_medium)
	# vE_Acs = reaction_equations.vE_Acs_Reaction()
	# vE_Acs.parameters={'kAcs_cat': 129000.0000000000,'KAcs_ACE': 0.0240000000}
	# model.append(vE_Acs)
	# vE_Acs_medium = reaction_equations.vE_Acs_medium_Reaction()
	# vE_Acs_medium.parameters={'KAcs_ACE': 0.0240000000,'kAcs_cat': 129000.0000000000,'rho': 564.0000000000}
	# model.append(vE_Acs_medium)
	# vE_Cs = reaction_equations.vE_Cs_Reaction()
	# vE_Cs.parameters={'KCs_OAA': 0.0017600000,'KCs_OAA_AcCoA': 0.0000886000,'KCs_AcCoA': 0.0298000000,'KCs_aKG': 0.1870000000,'kCs_cat': 5430000.0000000000}
	# model.append(vE_Cs)
	# vE_Icdh = reaction_equations.vE_Icdh_Reaction()
	# vE_Icdh.parameters={'kIcdh_cat': 34320000.0000000000,'KIcdh_ICIT': 0.0002010000,'KIcdh_PEP': 0.0548000000,'icdh_b': 0.0000006000,'LIcdh': 92.6000000000,'nIcdh': 2.0000000000}
	# model.append(vE_Icdh)
	# vE_akgdh = reaction_equations.vE_akgdh_Reaction()
	# vE_akgdh.parameters={'Kakgdh_NAD_m': 0.0552000000,'Kakgdh_SUC_I': 2.0700000000,'Kakgdh_NADH_I': 0.0179000000,'Kakgdh_CoA_m': 0.0031000000,'Kakgdh_aKG_m': 0.2370000000,'Kakgdh_Z': 4.1700000000,'kakgdh_cat': 702000.0000000000,'Kakgdh_aKG_I': 1.0000000000}
	# model.append(vE_akgdh)
	# vE_Sdh = reaction_equations.vE_Sdh_Reaction()
	# vE_Sdh.parameters={'kSdh2_cat': 1960000.0000000000,'KSdh_eq': 29.7000000000,'KSdh_SUC_m': 0.1690000000,'kSdh1_cat': 1960000.0000000000}
	# model.append(vE_Sdh)
	# vE_Fum = reaction_equations.vE_Fum_Reaction()
	# vE_Fum.parameters={'KFum_eq': 12.9000000000,'kFum2_cat': 3820000.0000000000,'KFum_FUM_m': 0.0965000000,'kFum1_cat': 3820000.0000000000}
	# model.append(vE_Fum)
	# vE_Mdh = reaction_equations.vE_Mdh_Reaction()
	# vE_Mdh.parameters={'KMdh_OAA_m': 0.2080000000,'KMdh_MAL_m': 0.1770000000,'kMdh2_cat': 23000000.0000000000,'kMdh1_cat': 23000000.0000000000,'KMdh_NADH_I': 0.0248000000,'KMdh_NAD_II': 0.9450000000,'KMdh_OAA_I': 0.3440000000,'KMdh_OAA_II': 0.0732000000,'KMdh_eq': 0.7500000000,'KMdh_NAD_I': 0.1140000000,'KMdh_NAD_m': 0.0522000000,'KMdh_NADH_m': 0.0313000000,'KMdh_MAL_I': 0.3460000000}
	# model.append(vE_Mdh)
	# vE_MaeB = reaction_equations.vE_MaeB_Reaction()
	# vE_MaeB.parameters={'KMaeB_AcCoA': 1.8200000000,'nMaeB': 1.9900000000,'LMaeB': 239000.0000000000,'KMaeB_cAMP': 4.5500000000,'kMaeB_cat': 18400000.0000000000,'KMaeB_MAL': 0.0014600000}
	# model.append(vE_MaeB)
	# vE_Pck = reaction_equations.vE_Pck_Reaction()
	# vE_Pck.parameters={'KPck_ATP_I': 0.0398000000,'KPck_OAA_I': 0.3470000000,'KPck_ADP_i': 0.0198000000,'KPck_ATP_i': 0.0401000000,'KPck_PEP': 0.0701000000,'KPck_PEP_i': 0.0600000000,'kPck_cat': 34500000.0000000000,'KPck_OAA': 0.5780000000}
	# model.append(vE_Pck)
	# vE_Ppc = reaction_equations.vE_Ppc_Reaction()
	# vE_Ppc.parameters={'nPpc': 5.0000000000,'LPpc': 5650000.0000000000,'KPpc_PEP': 0.0272000000,'kPpc_cat': 485100.0000000000,'KPpc_FBP': 0.1890000000}
	# model.append(vE_Ppc)
	# vE_Icl = reaction_equations.vE_Icl_Reaction()
	# vE_Icl.parameters={'nIcl': 4.0000000000,'KIcl_PEP': 0.0163000000,'KIcl_aKG': 0.8340000000,'kIcl_cat': 94200000.0000000000,'KIcl_ICIT': 0.0075500000,'LIcl': 191000.0000000000,'KIcl_3PG': 0.5220000000}
	# model.append(vE_Icl)
	# vE_Ms = reaction_equations.vE_Ms_Reaction()
	# vE_Ms.parameters={'kMs_cat': 17000000000.0000000000,'KMs_GOX': 1.1100000000,'KMs_AcCoA': 0.0462000000,'KMs_GOX_AcCoA': 0.0396000000}
	# model.append(vE_Ms)
	# vE_AceKki = reaction_equations.vE_AceKki_Reaction()
	# vE_AceKki.parameters={'nAceK': 5.0000000000,'LAceK': 283000000.0000000000,'KAceK_Icdh': 0.1970000000,'kAceKki_cat': 10000000000000000.0000000000,'KAceK_OAA': 0.0100000000}
	# model.append(vE_AceKki)
	# vE_AceKph = reaction_equations.vE_AceKph_Reaction()
	# vE_AceKph.parameters={'KAceK_IcdhP': 7.2600000000,'nAceKph': 5.0000000000,'KAceKph_OAA': 0.0100000000,'LAceKph': 70000.0000000000,'kAceKph_cat': 600.0000000000}
	# model.append(vE_AceKph)
	# vE_G6pdh = reaction_equations.vE_G6pdh_Reaction()
	# vE_G6pdh.parameters={'KG6pdh_NADPH_nadpinh': 0.0815000000,'KG6pdh_NADPH_g6pinh': 0.1930000000,'KG6pdh_NADP': 0.0038200000,'KG6pdh_G6P': 0.3560000000,'vG6pdh_max': 16600.0000000000}
	# model.append(vE_G6pdh)
	# vE_Pgl = reaction_equations.vE_Pgl_Reaction()
	# vE_Pgl.parameters={'KPgl_h2': 0.0000097000,'vPgl_max': 22800.0000000000,'KPgl_6PG_m': 10.0000000000,'KPgl_6PGL_m': 0.0230000000,'KPgl_eq': 42.7000000000,'KPgl_h1': 0.0043700000}
	# model.append(vE_Pgl)
	# vE_Edd = reaction_equations.vE_Edd_Reaction()
	# vE_Edd.parameters={'pK_Edd': 8.7400000000,'KEdd_eq': 1000.0000000000,'KEdd_KDPG_m': 2.0200000000,'KEdd_6PG_m': 0.1180000000,'pH_Edd_m': 7.5300000000,'vEdd_max': 515.0000000000}
	# model.append(vE_Edd)
	# vE_Eda = reaction_equations.vE_Eda_Reaction()
	# vE_Eda.parameters={'pK_Eda': 37.0000000000,'vEda_max': 667.0000000000,'KEda_PYR_m': 7.7000000000,'KEda_KDPG_m': 0.1500000000,'KEda_GAP_m': 1.1800000000,'pH_Eda_m': 10.4000000000,'KEda_eq': 0.5020000000}
	# model.append(vE_Eda)
	# vE_6Pgdh = reaction_equations.vE_6Pgdh_Reaction()
	# vE_6Pgdh.parameters={'K6Pgdh_NADP': 0.0208000000,'K6Pgdh_ATP_inh': 3.0100000000,'v6Pgdh_max': 248000.0000000000,'K6Pgdh_6PG': 0.1010000000,'K6Pgdh_NADPH_inh': 0.0421000000}
	# model.append(vE_6Pgdh)
	# vE_R5pi = reaction_equations.vE_R5pi_Reaction()
	# vE_R5pi.parameters={'vR5pi_max': 32100.0000000000,'KR5pi_eq': 0.4770000000}
	# model.append(vE_R5pi)
	# vE_Rpe = reaction_equations.vE_Rpe_Reaction()
	# vE_Rpe.parameters={'KRpe_eq': 1.4100000000,'vRpe_max': 12700.0000000000}
	# model.append(vE_Rpe)
	# vE_Tkt1 = reaction_equations.vE_Tkt1_Reaction()
	# vE_Tkt1.parameters={'KTkt1_eq': 1.2000000000,'vTkt1_max': 8870.0000000000}
	# model.append(vE_Tkt1)
	# vE_Tkt2 = reaction_equations.vE_Tkt2_Reaction()
	# vE_Tkt2.parameters={'vTkt2_max': 380000.0000000000,'KTkt2_eq': 9.9700000000}
	# model.append(vE_Tkt2)
	# vE_Tal = reaction_equations.vE_Tal_Reaction()
	# vE_Tal.parameters={'vTal_max': 71700.0000000000,'KTal_eq': 1.0500000000}
	# model.append(vE_Tal)
	# vE_Cya = reaction_equations.vE_Cya_Reaction()
	# vE_Cya.parameters={'vCya_max': 9.4500000000,'KCya_EIIAP': 0.0020800000}
	# model.append(vE_Cya)
	# vE_cAMPdegr = reaction_equations.vE_cAMPdegr_Reaction()
	# vE_cAMPdegr.parameters={'vcAMPdegr_max': 9.2100000000,'KcAMPdegr_cAMP': 0.0483000000}
	# model.append(vE_cAMPdegr)
	# OP_NADH = reaction_equations.OP_NADH_Reaction()
	# OP_NADH.parameters={'Kakgdh_NAD_m': 0.0552000000,'KGapdh_GAP': 0.1520000000,'KMdh_OAA_II': 0.0732000000,'KGapdh_PGP': 0.1340000000,'KPdh_NAD_m': 0.4010000000,'KPdh_CoA_m': 0.0048200000,'KMdh_OAA_m': 0.2080000000,'KMdh_MAL_m': 0.1770000000,'Kakgdh_SUC_I': 2.0700000000,'KGapdh_eq': 0.3000000000,'kMdh2_cat': 23000000.0000000000,'kMdh1_cat': 23000000.0000000000,'KPdh_i': 68.3000000000,'Kakgdh_aKG_m': 0.2370000000,'KMdh_eq': 0.7500000000,'KMdh_NAD_I': 0.1140000000,'KMdh_MAL_I': 0.3460000000,'Kakgdh_NADH_I': 0.0179000000,'kPdh_cat': 48300000.0000000000,'Kakgdh_Z': 4.1700000000,'KGapdh_NAD': 0.4500000000,'KMdh_NADH_m': 0.0313000000,'POratio': 3.4800000000,'kGapdh_cat': 5050000000.0000000000,'KMdh_NADH_I': 0.0248000000,'KMdh_NAD_II': 0.9450000000,'Kakgdh_CoA_m': 0.0031000000,'KPdh_AcCoA_m': 0.0080000000,'KMdh_OAA_I': 0.3440000000,'KPdh_NADH_m': 0.0471000000,'KGapdh_NADH': 0.0208000000,'kakgdh_cat': 702000.0000000000,'KMdh_NAD_m': 0.0522000000,'KPdh_PYR_m': 1.0000000000,'Kakgdh_aKG_I': 1.0000000000}
	# model.append(OP_NADH)
	# OP_FADH2 = reaction_equations.OP_FADH2_Reaction()
	# OP_FADH2.parameters={'POratio_prime': 1.4900000000,'KSdh_eq': 29.7000000000,'kSdh2_cat': 1960000.0000000000,'KSdh_SUC_m': 0.1690000000,'kSdh1_cat': 1960000.0000000000}
	# model.append(OP_FADH2)
	mu = reaction_equations.mu_Reaction()
	mu.parameters={'kATP': 0.00000000887}
	model.append(mu)
	vgrowth = reaction_equations.vgrowth_Reaction()
	vgrowth.parameters={'kATP': 0.00000000887}
	model.append(vgrowth)
	vG_glk = reaction_equations.vG_glk_Reaction()
	vG_glk.parameters = {'Kglk_Cra': 0.0000122000, 'KCraFBP': 29.4000000000, 'kATP': 0.0000088700,
						 'vglk_Cra_bound': 0.0188333333, 'vglk_Cra_unbound': 2.5000000000, 'nCraFBP': 2.0000000000,
						 'Cratotal': 0.3000000000}
	model.append(vG_glk)
	vG_pfkA = reaction_equations.vG_pfkA_Reaction()
	vG_pfkA.parameters = {'KpfkA_Cra': 0.0000098700, 'KCraFBP': 29.4000000000, 'vpfkA_Cra_bound': 0.0171666667,
						  'vpfkA_Cra_unbound': 0.6250000000, 'kATP': 0.0000088700, 'nCraFBP': 2.0000000000,
						  'Cratotal': 0.3000000000}
	model.append(vG_pfkA)
	vG_fbp = reaction_equations.vG_fbp_Reaction()
	vG_fbp.parameters = {'KCraFBP': 29.4000000000, 'kATP': 0.0000088700, 'nCraFBP': 2.0000000000,
						 'Cratotal': 0.3000000000, 'vfbp_Cra_bound': 0.0171666667, 'Kfbp_Cra': 0.0375000000,
						 'vfbp_Cra_unbound': 0.0000000000}
	model.append(vG_fbp)
	vG_fbaA = reaction_equations.vG_fbaA_Reaction()
	vG_fbaA.parameters = {'vfbaA_Cra_unbound': 0.0305000000, 'Crptotal': 11.5000000000, 'vfbaA_Crp_bound': 0.0188333333,
						  'vfbaA_Cra_bound': 0.0000000000, 'KCraFBP': 29.4000000000, 'KfbaA_Crp': 9.4800000000,
						  'kATP': 0.0000088700, 'nCraFBP': 2.0000000000, 'vfbaA_Crp_unbound': 0.0000000000,
						  'nCrpcAMP': 1.0000000000, 'Cratotal': 0.3000000000, 'KCrpcAMP': 420.0000000000,
						  'KfbaA_Cra': 3.2600000000}
	model.append(vG_fbaA)
	vG_gapA = reaction_equations.vG_gapA_Reaction()
	vG_gapA.parameters = {'Crptotal': 11.5000000000, 'vgapA_Crp_unbound': 0.0000000000, 'KCraFBP': 29.4000000000,
						  'kATP': 0.0000088700, 'nCraFBP': 2.0000000000, 'KgapA_Cra': 16.1000000000,
						  'nCrpcAMP': 1.0000000000, 'vgapA_Cra_unbound': 0.0266666667, 'KgapA_Crp': 47.5000000000,
						  'vgapA_Crp_bound': 0.0363333333, 'KCrpcAMP': 420.0000000000, 'Cratotal': 0.3000000000,
						  'vgapA_Cra_bound': 0.0000000000}
	model.append(vG_gapA)
	vG_pykF = reaction_equations.vG_pykF_Reaction()
	vG_pykF.parameters = {'KCraFBP': 29.4000000000, 'kATP': 0.0000088700, 'nCraFBP': 2.0000000000,
						  'KpykF_Cra': 0.0726000000, 'vpykF_Cra_bound': 0.0005200000, 'Cratotal': 0.3000000000,
						  'vpykF_Cra_unbound': 0.1366666667}
	model.append(vG_pykF)
	# vG_ppsA = reaction_equations.vG_ppsA_Reaction()
	# vG_ppsA.parameters={'KppsA_Cra': 0.0004880000,'nCraFBP': 2.0000000000,'kATP': 0.0000126000,'vppsA_Cra_bound': 0.0038600000,'Cratotal': 0.0003000000,'vppsA_Cra_unbound': 0.0000000000,'KCraFBP': 0.0294000000}
	# model.append(vG_ppsA)
	# vG_lpd = reaction_equations.vG_lpd_Reaction()
	# vG_lpd.parameters={'vlpd_PdhR_unbound': 0.0092900000,'PdhRtotal': 0.0000666000,'kATP': 0.0000126000,'KPdhRPYR': 0.0434000000,'Klpd_PdhR': 0.0000246000,'nPdhRPYR': 1.0000000000,'vlpd_PdhR_bound': 0.0000545000}
	# model.append(vG_lpd)
	# vG_acs = reaction_equations.vG_acs_Reaction()
	# vG_acs.parameters={'nCrpcAMP': 1.0000000000,'vacs_Crp_bound': 0.0026200000,'nacs': 2.3100000000,'kATP': 0.0000126000,'vacs_Crp_unbound': 0.0000000000,'KCrpcAMP': 0.4200000000,'Kacs_Crp': 0.0013600000,'Crptotal': 0.0115000000}
	# model.append(vG_acs)
	# vG_gltA = reaction_equations.vG_gltA_Reaction()
	# vG_gltA.parameters={'nCrpcAMP': 1.0000000000,'kATP': 0.0000126000,'ngltA': 1.0700000000,'KCrpcAMP': 0.4200000000,'Crptotal': 0.0115000000,'vgltA_Crp_bound': 0.0205000000,'KgltA_Crp': 0.0565000000,'vgltA_Crp_unbound': 0.0000000000}
	# model.append(vG_gltA)
	# vG_icdA = reaction_equations.vG_icdA_Reaction()
	# vG_icdA.parameters={'KicdA_Cra': 0.0000292000,'vicdA_Cra_unbound': 0.0005130000,'kATP': 0.0000126000,'nCraFBP': 2.0000000000,'Cratotal': 0.0003000000,'KCraFBP': 0.0294000000,'vicdA_Cra_bound': 0.0014200000}
	# model.append(vG_icdA)
	# vG_sucAB = reaction_equations.vG_sucAB_Reaction()
	# vG_sucAB.parameters={'nCrpcAMP': 1.0000000000,'nsucAB': 0.7400000000,'vsucAB_Crp_unbound': 0.0000000000,'kATP': 0.0000126000,'vsucAB_Crp_bound': 0.1740000000,'KCrpcAMP': 0.4200000000,'KsucAB_Crp': 0.3090000000,'Crptotal': 0.0115000000}
	# model.append(vG_sucAB)
	# vG_sdhCDAB = reaction_equations.vG_sdhCDAB_Reaction()
	# vG_sdhCDAB.parameters={'nCrpcAMP': 1.0000000000,'vsdhCDAB_Crp_unbound': 0.0000000000,'KsdhCDAB_Crp': 0.0856000000,'vsdhCDAB_Crp_bound': 0.0177000000,'kATP': 0.0000126000,'KCrpcAMP': 0.4200000000,'nsdhCDAB': 0.7400000000,'Crptotal': 0.0115000000}
	# model.append(vG_sdhCDAB)
	# vG_fumABC = reaction_equations.vG_fumABC_Reaction()
	# vG_fumABC.parameters={'nfumABC': 0.7400000000,'nCrpcAMP': 1.0000000000,'KfumABC_Crp': 0.1490000000,'kATP': 0.0000126000,'KCrpcAMP': 0.4200000000,'vfumABC_Crp_bound': 0.0376000000,'vfumABC_Crp_unbound': 0.0000000000,'Crptotal': 0.0115000000}
	# model.append(vG_fumABC)
	# vG_mdh = reaction_equations.vG_mdh_Reaction()
	# vG_mdh.parameters={'vmdh_Crp_unbound': 0.0000000000,'nCrpcAMP': 1.0000000000,'kATP': 0.0000126000,'vmdh_Crp_bound': 0.0733000000,'KCrpcAMP': 0.4200000000,'Kmdh_Crp': 0.2140000000,'Crptotal': 0.0115000000}
	# model.append(vG_mdh)
	# vG_maeB = reaction_equations.vG_maeB_Reaction()
	# vG_maeB.parameters={'SS_MaeB': 0.0075000000,'kATP': 0.0000126000}
	# model.append(vG_maeB)
	# vG_pckA = reaction_equations.vG_pckA_Reaction()
	# vG_pckA.parameters={'KpckA_Cra': 0.0001520000,'nCraFBP': 2.0000000000,'kATP': 0.0000126000,'vpckA_Cra_unbound': 0.0000000000,'vpckA_Cra_bound': 0.0099500000,'Cratotal': 0.0003000000,'KCraFBP': 0.0294000000}
	# model.append(vG_pckA)
	# vG_ppc = reaction_equations.vG_ppc_Reaction()
	# vG_ppc.parameters={'SS_Ppc': 0.0074000000,'kATP': 0.0000126000}
	# model.append(vG_ppc)
	# vG_aceA = reaction_equations.vG_aceA_Reaction()
	# vG_aceA.parameters={'IclRtotal': 0.0000830000,'Kmu': 0.6000000000,'KaceBAK_Crp': 0.5260000000,'vaceBAK_Cra_unbound': 0.0000127000,'n_mu': 4.0000000000,'KaceBAK_GOX': 0.0022500000,'KaceBAK_Cra': 0.0005400000,'kATP': 0.0000126000,'kaceBAK_cat_IclR': 0.3200000000,'vaceBAK_Crp_bound': 0.0000015800,'vaceBAK_Cra_bound': 0.0037000000,'nCraFBP': 2.0000000000,'KaceBAK_PYRprime': 0.0027200000,'KCrpcAMP': 0.4200000000,'Cratotal': 0.0003000000,'KCraFBP': 0.0294000000,'nCrpcAMP': 1.0000000000,'vaceBAK_Crp_unbound': 0.0002510000,'KaceBAK_DNA': 0.0000009510,'LaceBAK': 412.0000000000,'KaceBAK_PYR': 2.0100000000,'Crptotal': 0.0115000000}
	# model.append(vG_aceA)
	# vG_aceB = reaction_equations.vG_aceB_Reaction()
	# vG_aceB.parameters={'IclRtotal': 0.0000830000,'Kmu': 0.6000000000,'KaceBAK_Crp': 0.5260000000,'vaceBAK_Cra_unbound': 0.0000127000,'n_mu': 4.0000000000,'Factor_aceB': 0.3090000000,'KaceBAK_GOX': 0.0022500000,'KaceBAK_Cra': 0.0005400000,'kATP': 0.0000126000,'kaceBAK_cat_IclR': 0.3200000000,'vaceBAK_Crp_bound': 0.0000015800,'vaceBAK_Cra_bound': 0.0037000000,'nCraFBP': 2.0000000000,'KaceBAK_PYRprime': 0.0027200000,'KCrpcAMP': 0.4200000000,'Cratotal': 0.0003000000,'KCraFBP': 0.0294000000,'nCrpcAMP': 1.0000000000,'vaceBAK_Crp_unbound': 0.0002510000,'KaceBAK_DNA': 0.0000009510,'LaceBAK': 412.0000000000,'KaceBAK_PYR': 2.0100000000,'Crptotal': 0.0115000000}
	# model.append(vG_aceB)
	# vG_aceK = reaction_equations.vG_aceK_Reaction()
	# vG_aceK.parameters={'IclRtotal': 0.0000830000,'Kmu': 0.6000000000,'KaceBAK_Crp': 0.5260000000,'vaceBAK_Cra_unbound': 0.0000127000,'n_mu': 4.0000000000,'KaceBAK_GOX': 0.0022500000,'KaceBAK_Cra': 0.0005400000,'kATP': 0.0000126000,'kaceBAK_cat_IclR': 0.3200000000,'Factor_aceK': 0.0172000000,'vaceBAK_Crp_bound': 0.0000015800,'vaceBAK_Cra_bound': 0.0037000000,'nCraFBP': 2.0000000000,'KaceBAK_PYRprime': 0.0027200000,'KCrpcAMP': 0.4200000000,'Cratotal': 0.0003000000,'KCraFBP': 0.0294000000,'nCrpcAMP': 1.0000000000,'vaceBAK_Crp_unbound': 0.0002510000,'KaceBAK_DNA': 0.0000009510,'LaceBAK': 412.0000000000,'KaceBAK_PYR': 2.0100000000,'Crptotal': 0.0115000000}
	# model.append(vG_aceK)
	# vD_X = reaction_equations.vD_X_Reaction()
	# vD_X.parameters={'D': 0.0000000000}
	# model.append(vD_X)
	# vD_GLCfeed = reaction_equations.vD_GLCfeed_Reaction()
	# vD_GLCfeed.parameters={'D': 0.0000000000}
	# model.append(vD_GLCfeed)
	# vD_GLCex = reaction_equations.vD_GLCex_Reaction()
	# vD_GLCex.parameters={'D': 0.0000000000}
	# model.append(vD_GLCex)
	# vD_GLC = reaction_equations.vD_GLC_Reaction()
	# vD_GLC.parameters={'kATP': 0.0000126000}
	# model.append(vD_GLC)
	# vD_G6P = reaction_equations.vD_G6P_Reaction()
	# vD_G6P.parameters={'kATP': 0.0000126000}
	# model.append(vD_G6P)
	# vD_F6P = reaction_equations.vD_F6P_Reaction()
	# vD_F6P.parameters={'kATP': 0.0000126000}
	# model.append(vD_F6P)
	# vD_FBP = reaction_equations.vD_FBP_Reaction()
	# vD_FBP.parameters={'kATP': 0.0000126000}
	# model.append(vD_FBP)
	# vD_GAP = reaction_equations.vD_GAP_Reaction()
	# vD_GAP.parameters={'kATP': 0.0000126000}
	# model.append(vD_GAP)
	# vD_PEP = reaction_equations.vD_PEP_Reaction()
	# vD_PEP.parameters={'kATP': 0.0000126000}
	# model.append(vD_PEP)
	# vD_PYR = reaction_equations.vD_PYR_Reaction()
	# vD_PYR.parameters={'kATP': 0.0000126000}
	# model.append(vD_PYR)
	# vD_AcCoA = reaction_equations.vD_AcCoA_Reaction()
	# vD_AcCoA.parameters={'kATP': 0.0000126000}
	# model.append(vD_AcCoA)
	# vD_AcP = reaction_equations.vD_AcP_Reaction()
	# vD_AcP.parameters={'kATP': 0.0000126000}
	# model.append(vD_AcP)
	# vD_ACEex = reaction_equations.vD_ACEex_Reaction()
	# vD_ACEex.parameters={'D': 0.0000000000}
	# model.append(vD_ACEex)
	# vD_ICIT = reaction_equations.vD_ICIT_Reaction()
	# vD_ICIT.parameters={'kATP': 0.0000126000}
	# model.append(vD_ICIT)
	# vD_aKG = reaction_equations.vD_aKG_Reaction()
	# vD_aKG.parameters={'kATP': 0.0000126000}
	# model.append(vD_aKG)
	# vD_SUC = reaction_equations.vD_SUC_Reaction()
	# vD_SUC.parameters={'kATP': 0.0000126000}
	# model.append(vD_SUC)
	# vD_FUM = reaction_equations.vD_FUM_Reaction()
	# vD_FUM.parameters={'kATP': 0.0000126000}
	# model.append(vD_FUM)
	# vD_MAL = reaction_equations.vD_MAL_Reaction()
	# vD_MAL.parameters={'kATP': 0.0000126000}
	# model.append(vD_MAL)
	# vD_OAA = reaction_equations.vD_OAA_Reaction()
	# vD_OAA.parameters={'kATP': 0.0000126000}
	# model.append(vD_OAA)
	# vD_GOX = reaction_equations.vD_GOX_Reaction()
	# vD_GOX.parameters={'kATP': 0.0000126000}
	# model.append(vD_GOX)
	# vD_6PGL = reaction_equations.vD_6PGL_Reaction()
	# vD_6PGL.parameters={'kATP': 0.0000126000}
	# model.append(vD_6PGL)
	# vD_6PG = reaction_equations.vD_6PG_Reaction()
	# vD_6PG.parameters={'kATP': 0.0000126000}
	# model.append(vD_6PG)
	# vD_KDPG = reaction_equations.vD_KDPG_Reaction()
	# vD_KDPG.parameters={'kATP': 0.0000126000}
	# model.append(vD_KDPG)
	# vD_RU5P = reaction_equations.vD_RU5P_Reaction()
	# vD_RU5P.parameters={'kATP': 0.0000126000}
	# model.append(vD_RU5P)
	# vD_R5P = reaction_equations.vD_R5P_Reaction()
	# vD_R5P.parameters={'kATP': 0.0000126000}
	# model.append(vD_R5P)
	# vD_X5P = reaction_equations.vD_X5P_Reaction()
	# vD_X5P.parameters={'kATP': 0.0000126000}
	# model.append(vD_X5P)
	# vD_S7P = reaction_equations.vD_S7P_Reaction()
	# vD_S7P.parameters={'kATP': 0.0000126000}
	# model.append(vD_S7P)
	# vD_E4P = reaction_equations.vD_E4P_Reaction()
	# vD_E4P.parameters={'kATP': 0.0000126000}
	# model.append(vD_E4P)
	# vD_cAMP = reaction_equations.vD_cAMP_Reaction()
	# vD_cAMP.parameters={'kATP': 0.0000126000}
	# model.append(vD_cAMP)
	vD_Glk = reaction_equations.vD_Glk_Reaction()
	vD_Glk.parameters={'kdegr': 0.03300000000,'kATP': 0.00000000887}
	model.append(vD_Glk)

	vD_G6P = reaction_equations.vD_G6P_Reaction()
	vD_G6P.parameters = {'kATP': 0.0000088700}
	model.append(vD_G6P)
	vD_F6P = reaction_equations.vD_F6P_Reaction()
	vD_F6P.parameters = {'kATP': 0.0000088700}
	model.append(vD_F6P)
	vD_FBP = reaction_equations.vD_FBP_Reaction()
	vD_FBP.parameters = {'kATP': 0.0000088700}
	model.append(vD_FBP)
	vD_GAP = reaction_equations.vD_GAP_Reaction()
	vD_GAP.parameters = {'kATP': 0.0000088700}
	model.append(vD_GAP)
	vD_PEP = reaction_equations.vD_PEP_Reaction()
	vD_PEP.parameters = {'kATP': 0.0000088700}
	model.append(vD_PEP)
	vD_PYR = reaction_equations.vD_PYR_Reaction()
	vD_PYR.parameters = {'kATP': 0.0000088700}
	model.append(vD_PYR)
	# vD_Pps = reaction_equations.vD_Pps_Reaction()
	# vD_Pps.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Pps)
	# vD_Pdh = reaction_equations.vD_Pdh_Reaction()
	# vD_Pdh.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Pdh)
	# vD_Acs = reaction_equations.vD_Acs_Reaction()
	# vD_Acs.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Acs)
	# vD_Cs = reaction_equations.vD_Cs_Reaction()
	# vD_Cs.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Cs)
	# vD_Icdh = reaction_equations.vD_Icdh_Reaction()
	# vD_Icdh.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Icdh)
	# vD_IcdhP = reaction_equations.vD_IcdhP_Reaction()
	# vD_IcdhP.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_IcdhP)
	# vD_akgdh = reaction_equations.vD_akgdh_Reaction()
	# vD_akgdh.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_akgdh)
	# vD_Sdh = reaction_equations.vD_Sdh_Reaction()
	# vD_Sdh.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Sdh)
	# vD_Fum = reaction_equations.vD_Fum_Reaction()
	# vD_Fum.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Fum)
	# vD_Mdh = reaction_equations.vD_Mdh_Reaction()
	# vD_Mdh.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Mdh)
	# vD_MaeB = reaction_equations.vD_MaeB_Reaction()
	# vD_MaeB.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_MaeB)
	# vD_Pck = reaction_equations.vD_Pck_Reaction()
	# vD_Pck.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Pck)
	# vD_Ppc = reaction_equations.vD_Ppc_Reaction()
	# vD_Ppc.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Ppc)
	# vD_Icl = reaction_equations.vD_Icl_Reaction()
	# vD_Icl.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Icl)
	# vD_Ms = reaction_equations.vD_Ms_Reaction()
	# vD_Ms.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_Ms)
	# vD_AceK = reaction_equations.vD_AceK_Reaction()
	# vD_AceK.parameters={'kdegr': 0.3300000000,'kATP': 0.0000126000}
	# model.append(vD_AceK)
	# vBM_G6P = reaction_equations.vBM_G6P_Reaction()
	# vBM_G6P.parameters = {'kBM_GLC_G6P': 0.8150000000}
	# model.append(vBM_G6P)
	# vBM_F6P = reaction_equations.vBM_F6P_Reaction()
	# vBM_F6P.parameters = {'kBM_GLC_F6P': 21.1666666667}
	# model.append(vBM_F6P)
	# vBM_GAP = reaction_equations.vBM_GAP_Reaction()
	# vBM_GAP.parameters = {'kBM_GLC_GAP': 1.7000000000}
	# model.append(vBM_GAP)
	# vBM_PEP = reaction_equations.vBM_PEP_Reaction()
	# vBM_PEP.parameters = {'kBM_GLC_PEP': 15.9500000000}
	# model.append(vBM_PEP)
	# vBM_PYR = reaction_equations.vBM_PYR_Reaction()
	# vBM_PYR.parameters = {'kBM_GLC_PYR': 14.1833333333}
	# model.append(vBM_PYR)
	# vBM_AcCoA = reaction_equations.vBM_AcCoA_Reaction()
	# vBM_AcCoA.parameters={'kBM_GLC_AcCoA': 1659.0000000000}
	# model.append(vBM_AcCoA)
	# vBM_aKG = reaction_equations.vBM_aKG_Reaction()
	# vBM_aKG.parameters={'kBM_GLC_aKG': 300.0000000000}
	# model.append(vBM_aKG)
	# vBM_SUC = reaction_equations.vBM_SUC_Reaction()
	# vBM_SUC.parameters={'kBM_GLC_SUC': 190.0000000000}
	# model.append(vBM_SUC)
	# vBM_FUM = reaction_equations.vBM_FUM_Reaction()
	# vBM_FUM.parameters={'kBM_GLC_FUM': 3470.0000000000}
	# model.append(vBM_FUM)
	# vBM_OAA = reaction_equations.vBM_OAA_Reaction()
	# vBM_OAA.parameters={'kBM_GLC_OAA': 23100.0000000000}
	# model.append(vBM_OAA)
	# vBM_R5P = reaction_equations.vBM_R5P_Reaction()
	# vBM_R5P.parameters={'kBM_GLC_R5P': 308.0000000000}
	# model.append(vBM_R5P)
	# vBM_E4P = reaction_equations.vBM_E4P_Reaction()
	# vBM_E4P.parameters={'kBM_GLC_E4P': 1510.0000000000}
	# model.append(vBM_E4P)
	model.set_time(0, 300, 1000)

	model.species = {"X": 12.000000,
		 "GLCex": 22500.000000,
		 "GLC": 2845.000000,
		 "G6P": 3479.000000,
		 "F6P" : 600,
    "FBP" : 15000,
	"BPG": 65.425903,
	"PGA3": 696.540882,
	"PGA2": 378.442756,
    "GAP" : 218.200,
    "PEP" : 402.800,
    "PYR" : 1765.000000,
    "AcCoA" : 0.611700,
    "AcP" : 1.103000,
    "ACEex" : 0.000000,
    "ICIT" : 0.109900,
    "aKG" : 0.441800,
    "SUC" : 0.573400,
    "FUM" : 0.119900,
    "MAL" : 1.432000,
    "OAA" : 0.650700,
    "GOX" : 0.080510,
    "sixPGL" : 1607.000,
    "sixPG" : 3457.000,
    "KDPG" : 1.041000,
    "RU5P" : 0.158600,
    "R5P" : 0.397900,
    "X5P" : 0.187800,
    "S7P" : 0.066270,
    "E4P" : 0.097020,
    "cAMP" : 0.019630,
    "EIIAP" : 0.014450,

    "Glk" : 0.003066,
    "Pfk" : 0.001432,
    "Fbp" : 0.000297,

    "Fba" : 0.023420,
    "Gapdh" : 0.040280,
    "Pyk" : 0.002879,
    "Pps" : 0.003863,
    "Pdh" : 0.000530,
    "Acs" : 0.000169,
    "Cs" : 0.002777,
    "Icdh" : 0.172100,
    "IcdhP" : 0.119500,
    "akgdh" : 0.002612,
    "Sdh" : 0.059540,
    "Fum" : 0.002105,
    "Mdh" : 0.007350,
    "MaeB" : 0.000187,
    "Pck" : 0.000673,
    "Ppc" : 0.000003,
    "Icl" : 0.012040,
    "Ms7" : 0.003618,
    "AceK" : 0.000029,

    "ATP" : 9600.000,
    "ADP" : 560.000,

    "AMP" : 0.280000,
    "NAD" : 2.600000,
    "NADH" : 0.083000,
    "NADP" : 0.002100,
    "NADPH" : 0.120000,
    "CoA" : 1.400000,
	}

	model.setup_model()
	model.run_model()


	model.plot_substrate("GLC",plot=True)
	model.plot_substrate("G6P", plot=True)
	model.plot_substrate("F6P", plot=True)
	model.plot_substrate("FBP", plot=True)
	model.plot_substrate("GAP", plot=True)
	model.plot_substrate("DAP", plot=True)
	model.plot_substrate("BPG", plot=True)
	model.plot_substrate("PGA3", plot=True)
	model.plot_substrate("PGA2", plot=True)
	model.plot_substrate("PEP", plot=True)
	model.plot_substrate("PYR", plot=True)
	# model.plot_substrate("AcCoa", plot=True)
	# model.plot_substrate("ICIT", plot=True)
	# model.plot_substrate("aKG", plot=True)
	# model.plot_substrate("SUC", plot=True)
	# model.plot_substrate("FUM", plot=True)
	# model.plot_substrate("MAL", plot=True)
	# model.plot_substrate("OAA", plot=True)