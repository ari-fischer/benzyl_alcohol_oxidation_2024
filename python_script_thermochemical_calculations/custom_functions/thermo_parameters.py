def collect_E_empty(df):
#takes in a df, and then finds the important energy information contained in the paths
    import numpy as np
    import os
    Es = np.array([])
    G_cds = np.array([])
    for i in range(df.shape[0]):
        lp = df['path_l_geom'][i]
        print(lp)
        if str(lp) == 'nan':
            E = 0
            G = 0
        else:
            #import energy calculations
            if os.path.exists(lp+'/Energy.txt'):
                with open(lp+'/Energy.txt') as f:
                    lines = f.readlines() 
                if np.size(lines) > 0:
                    E = float(lines[0].split()[3])
                else:
                    E = 0
                #import cavitation, dispersion, solvation   
                with open(lp+'/G_cds.txt') as f:
                    lines = f.readlines() 
                if np.size(lines) > 0:
                    #print(lines[0].split())
                    G = float(lines[0].split()[2])
                else:
                    G = 0
            else:
                E = 0
                G = 0
        Es = np.append(Es,E)
        G_cds = np.append(G_cds,G)
    
    #add to the dataframe
    #df.insert(2,'E',Es*2625.5) #kJ/mol
    df['G_cds'] = G_cds*4.1844

    #df.insert(3,'G_cds',G_cds*4.1844) #kcal/mol -> kJ    

    #if 'E' in df.columns:
    df['E'] = Es*2625.5
    #    df.insert(3,'E',Es*2625.5) #Ha -> kJ 
    return df


def collect_freq_empty(df,T_C):

#takes in dataframe from CSV with all the important info
#updated to export # of imaginary and the value of the first one

#import information from frequency calculations
# This block calculates the thermodynamic information based on vibrational freuencies as harmonic oscillators

#July 2024- imaginary export frequency
    import numpy as np
    import pandas as pd
    
    [N_A,k_B,T,h,kJmol_Ha] = [6.0221408e23,1.380649e-23, 273+T_C,6.626070E-34,2625.5] #[ -- , m2 kg s-2 K-1, K,  # m2 kg / s]
    R_ig = 8.3144E-3
    R = k_B*N_A # ideal gas constant: J / mol / K


    h_ev = 4.135668E-15 #planks constant in ev*s
    k_ev = 0.0000861733 #boltzmann in ev/K

    #initial vector to build S H and ZPEs
    S_vec = np.array([])
    H_vec = np.array([])
    ZPE_vec = np.array([])
    H_rot_vec = np.array([])
    H_trans_vec = np.array([])
    S_rot_ig_vec = np.array([])
    S_trans_vec = np.array([])
    S_elec_vec = np.array([])
    n_ims = np.array([])
    vibs_list = np.array([])
    m_red = np.array([])
    f_ims = np.array([])

    for i in range(df.shape[0]):
        
        lp = df['path_l_freq'][i]
        print(lp)
        print(df['Calc'][i])
        if str(lp) == 'nan':
            S_vec = np.append(S_vec,0) 
            #Data['H_vib'][i] = np.sum(np.sum(S_damp))
            H_vec = np.append(H_vec,0) 
            ZPE_vec = np.append(ZPE_vec,0)
            
            H_rot_vec = np.append(H_rot_vec,0)
            H_trans_vec = np.append(H_trans_vec,0)
            S_rot_ig_vec = np.append(S_rot_ig_vec,0)
            S_trans_vec = np.append(S_trans_vec,0)
            S_elec_vec = np.append(S_elec_vec,0)
            n_ims = np.append(n_ims,[0]) 
        else:
            
            with open(lp+'/Frequencies.txt') as f:
                lines = f.readlines() 
            if np.size(lines) > 0: #checks that the freq completed and gave answer output

                vibs = pd.read_csv(lp+'/Frequencies.txt',skiprows=0,header=None,delim_whitespace=True)
                vib_0 = vibs[0][0]
                #moments of iniertial
                #flag if there are multiple 
                if df['Calc'][i][0:2]=='TS':
                    #for TS excluding one of the modes corresponding to the direction along the rxn coordiante. TS needs to be in the folder name
                    n_im = np.sum(vibs[1]<0)
                else:
                    n_im = np.sum(vibs[0]<0)
                n_ims = np.append(n_ims,n_im)
                f_ims = np.append(f_ims,vib_0)
                
                #get moment of inertia    
                I_iner = pd.read_csv(lp+'/inertia.txt',skiprows=0,header=None,delim_whitespace=True).to_numpy()
                I_iner = I_iner[0] #amu / bohr

                #get molecular mass
                M_mass = pd.read_csv(lp+'/M_mass.txt',skiprows=0,header=None,delim_whitespace=True).to_numpy()
                m = M_mass[0][2]*1.67377E-27 #amu to kg

                #get symmetry number:
                sym = pd.read_csv(lp+'/Sym.txt',skiprows=0,header=None,delim_whitespace=True).to_numpy()
                sigma = sym[0][4] #amu / bohr

                print(sigma)

                linear = df['Smiles'][i] in list(['[O][O-]','[O][O]','[OH]','[OH-]'])
                #find the one principle moment of inertia for a linear molecule
                if linear:
                    I_iner = np.max(I_iner)

                #apply threshold from grimme and others for the effect of low vib modes on entropy and enthalpy
                #vibs[vibs<100]=100 #in cm01
                print(vibs)
                #convert to hz
                v_hz = vibs.to_numpy()*29979.2458*1e6#cm-1 -> Mhz -> hz
                #damppend vibrations as per Grimme... damping and cutoff at 100?
                alpha = 4
                w_vib = 1/(1+(100/vibs)**alpha)
                #calculate vibrational modes
                S_vibs =R*((h_ev*v_hz/(k_ev*T*(np.exp(h_ev*v_hz/k_ev/T)-1))-np.log(1-np.exp(-h_ev*v_hz/k_ev/T))))
                print(S_vibs)
                mu = h/8/v_hz/np.pi**2 # ev * s**2
                #amu*Bohr^2... convert to kJ/m2 then to ev

                I_ave = (np.mean(I_iner))*(5.29177210903e-11**2)*1.66053907e-27 #amu/bohr...
                mu_p = (mu*I_ave/(mu+I_ave))*6.24150907e18 # m2 kg -> ev...
                # dampened vibrationals/rottions for low wavenumber modes
                S_rot_damp = R*(1/2+np.log((8*np.pi**3*mu_p*k_ev*T/h_ev**2)**(1/2)))
                S_damp =w_vib*S_vibs+(1-w_vib)*S_rot_damp
                print(S_damp)
                #S_damp = S_vibs
                #calculate ZPE and vib enthalpy
                ZPE = 1/2*h_ev*v_hz
                H_vib = h_ev*v_hz*np.exp(-h_ev*v_hz/k_ev/T)/(1-np.exp(-h_ev*v_hz/k_ev/T))
                #populate table with S H and ZPE values:
                
                
                if df['Calc'][i][0:2]=='TS':
                    #for TS excluding one of the modes corresponding to the direction along the rxn coordiante. TS needs to be in the folder name
                    S_vec = np.append(S_vec,np.sum(np.sum(S_damp[1:][:]))+np.sum(S_damp[0][1:]))#
                    H_vec = np.append(H_vec,np.sum(np.sum(H_vib[1:][:]))+np.sum(H_vib[0][1:]))#
                    ZPE_vec = np.append(ZPE_vec,np.sum(np.sum(ZPE[1:][:]))+np.sum(ZPE[0][1:]))#
                    
                else:
                    #for minima keep all modes
                    S_vec = np.append(S_vec,np.sum(np.sum(S_damp))) 
                    #Data['H_vib'][i] = np.sum(np.sum(S_damp))
                    H_vec = np.append(H_vec,np.sum(np.sum(H_vib))) 
                    ZPE_vec = np.append(ZPE_vec,np.sum(np.sum(ZPE))) 
                
                # H and S rot and trans
                df_H = pd.read_csv(lp + '/H_out'  +'.txt',skiprows=0,header=None,delim_whitespace=True)*4.184
                # trans, rot, vib, tot
                df_S = pd.read_csv(lp + '/S_out'  +'.txt',skiprows=0,header=None,delim_whitespace=True)*4.184
                
                #populate the table
                #https://cccbdb.nist.gov/thermox.asp



                if linear:
                    B = h/8/np.pi/np.pi/(I_iner*(5.29177210903e-11**2)*1.66053907e-27)
                    H_rot = R_ig*T # kJ/mol 
                    #for polyatomic gas (diff for H)
                    #joules
                    S_rot_ig = R_ig*1000*(np.log(k_B*T/h/sigma/B)+1)
                    
                else: #polyatomic
                    H_rot = 3/2*R_ig*T # kJ/mol 
                    
                    
                    A = h/8/np.pi/np.pi/(I_iner[0]*(5.29177210903e-11**2)*1.66053907e-27)
                    B = h/8/np.pi/np.pi/(I_iner[1]*(5.29177210903e-11**2)*1.66053907e-27)
                    C = h/8/np.pi/np.pi/(I_iner[2]*(5.29177210903e-11**2)*1.66053907e-27)
                    H_rot = 3/2*R_ig*T# kJ/mol   
                    #joules 
                    S_rot_ig = R_ig*1000*(3/2*np.log(k_B*T/h)-1/2*np.log(A*B*C/np.pi)-np.log(sigma)+3/2)
                #joules
                S_trans = R_ig*1000*(3/2*np.log(2*np.pi*m/h/h)+5/2*np.log(k_B*T)-np.log(1E5)+5/2) # kJ/mol/K  
                #kilojoules
                H_trans = 5/2*R_ig*T #df_H.to_numpy()[0] # kJ/mol
                

                with open(df['path_l_geom'][i]+'/out_head.txt') as f:
                    lines = f.readlines()
                M = int(lines[-1][-2])        
            
                if df['Calc'][i]=='OH':
                    S_elec = R_ig*1000*np.log(4)
                else:
                    S_elec = R_ig*1000*np.log(M)

                
                H_rot_vec = np.append(H_rot_vec,H_rot)
                H_trans_vec = np.append(H_trans_vec,H_trans)
                S_rot_ig_vec = np.append(S_rot_ig_vec,S_rot_ig)
                S_trans_vec = np.append(S_trans_vec,S_trans)
                S_elec_vec = np.append(S_elec_vec,S_elec)
                
            else: #populates with zero if there is no computed frequencies
                S_vec = np.append(S_vec,0) 
                #Data['H_vib'][i] = np.sum(np.sum(S_damp))
                H_vec = np.append(H_vec,0) 
                ZPE_vec = np.append(ZPE_vec,0)
                
                H_rot_vec = np.append(H_rot_vec,0)
                H_trans_vec = np.append(H_trans_vec,0)
                S_rot_ig_vec = np.append(S_rot_ig_vec,0)
                S_trans_vec = np.append(S_trans_vec,0)
                S_elec_vec = np.append(S_elec_vec,0)
                n_ims = np.append(n_ims,[0]) 
                f_ims = np.append(f_ims,[0]) 

            #now read masses file . Added 28 Feb 2024     
            with open(lp+'/masses.txt') as f:
                lines = f.readlines() 
                if np.size(lines) > 0:
                    masses = pd.read_csv(lp+'/masses.txt',skiprows=0,header=None,delim_whitespace=True)
                    mass_0 = masses[0][0]
                else:
                    mass_0 = [0]
            
        m_red = np.append(m_red,mass_0)
        vibs_list = np.append(vibs_list,vib_0)
            


    df['vib_0'] = vibs_list   
    df['n_imag'] = n_ims
    df['f_imag'] = f_ims
    df['m_red'] = m_red
    df['S_vib'] = S_vec
    df['H_vib'] = H_vec*96.4869
    df['ZPE']= ZPE_vec*96.4869 #kJ/mol


    df['H_rot'] = H_rot_vec #kcal/mol -> kJ/mol (occured at import)
    df['H_trans'] = H_trans_vec #kcal/mol -> kJ/mol
    df['S_rot'] = S_rot_ig_vec #kcal/mol -> J/mol/K
    #modified 29 feb 24-- add in the change in standard state for the ln(p) in gas formalism
    # convert from 1 atm to 23.45309 atm (for 1 molar reference state)
    df['S_trans'] = S_trans_vec-np.log(23.4530933)*R_ig*1000 #kcal/mol -> J/mol/K 
    df['S_elec'] = S_elec_vec #kcal/mol -> J/mol/K




    return df


def species_list(s):
    s = [item for sublist in s for item in sublist]

    species = []
    for item in s:
        if '.' in item:
            species.extend(item.split('.'))
        else:
            species.append(item)
        
    species = list(set(species))
    #for item in species:
    #    print(item)
    return species

def rxn_strings(rcts,prod):
    s_rxns = list([])
    for i in range(len(rcts)):
        if len(rcts[i])==2:
            s_rcts = rcts[i][0]+'+'+rcts[i][1]
        elif len(rcts[i])==3:
            s_rcts = rcts[i][0]+'+'+rcts[i][1]+'+'+rcts[i][2]
        else:
            s_rcts = rcts[i][0]
        
        if len(prod[i])==2:
            s_prd = prod[i][0]+'+'+prod[i][1]
        elif len(prod[i])==3:
            s_prd = prod[i][0]+'+'+prod[i][1]+'+'+prod[i][2]
        elif len(prod[i])==4:
            s_prd = prod[i][0]+'+'+prod[i][1]+'+'+prod[i][2]+'+'+prod[i][3]    
        else:
            s_prd = prod[i][0]
        
        #print(s_rcts+'->'+s_prd)
        s_rxns.append(s_rcts+'->'+s_prd)
    return s_rxns