#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 13:29:22 2018

@author: boss 
"""
 
""" This is just a change to the discription to test 
pushing a change in the code.  
 """

import cantera
import matplotlib.pyplot as plt 
import numpy as np
import math
import os
import datetime
import time
#import pickle
#import cProfile
#from statistics import mode
cantera.suppress_thermo_warnings()

T =np.linspace(600,1000,2)
P =np.linspace(1,30,2)
Phi = np.linspace(0.1,2,2)
CH4 = np.linspace(0.001,0.01,2)
rxns = np.linspace(0,324,325) #N + NO <=> N2 + O; N + O2 <=> NO + O; N + OH <=> NO + H
#(0,324,325)
rxns = np.linspace(0,324,325) 
interested_rxns = [177, 178, 179, 239]
SpecificSpecies = ['H2O', 'CH4', 'CO', 'CO2', 'NO']
endtime = 5
cut_time = 0 #in case a reading is too early 


def reac_sens():
    """ Function which saves simulation history.
    
    This function saves the different parameters through each time step
    per simulation. Time, Temperature, Pressure, Mole Fractions, 
    Sensitivities. It reduces the amount of information through the
    calling of the 'reduce_resolution' function. Lastly, it seperates all
    information into different lists to facilitate other calculations.
    
    Dummy Variables
    ---------------
    row : list
        List which collects all information from current time step.
    rows1 : list
        List which appends all information from 'row' to pass on to
        'reduce resolution' function
    rows : list
        List which contains the reduced amount of points returned from
        the 'reduce_resolution' function
    
    Appends
    --------
    t_T_P_AMol : list
        List which contains time, Temperature, Pressure, all species mole
        fractions per simulation.
    t_SMol : list
        List which contains time, and 'SpecificSpecies' mole fractions
        per simulation.
    t_AllSpeciesSens : list
        List which contains time, and all sensitivities of species with
        respect to reactions available in GRI-Mech 3.0, an optimized 
        mechanism designed to model natural gas combustion, including 
        NO formation and reburn chemistry. List is unique per simulation.
    All_time_Sens : list
        List which saves time and all species sensitivities with respect
        to reactions available in GRI-Mech 3.0 mechanism file. This list
        is a list of lists, in which each list corresponds to a different
        simulation or set of initial conditions.
    All_tTP_AMo : list
        This is another list which saves time, Temperature, Pressure, and
        all species mole fractions. This list also is a list of lists in
        which each individual list corresponds to a different simulation
        or set of initial conditions.
    
    """
    
    for i in rxns:
        reac.add_sensitivity_reaction(i)
        
    row = [None]*6
    rows1=[]
    rows=[]
    
    while sim.time < endtime and sim.time >= cut_time:
        row[0] = sim.step()          #save time
        row[1:3] = reac.thermo.TP     #save T, P
        row[3] = gas1[gas1.species_names].X # save mole fraction value of all species
        row[4] = gas1[SpecificSpecies].X   #save molesfrac value of specific species
        ssensitivities = sim.sensitivities() #sensitivities
        row[5] = ssensitivities[2:55] #first two rows are mass and enthalpy or temperature, not sensitivities
        rows1.append([x for x in row])
    rows = reduce_resolution(rows1,500) #isnt rows 1 always going to have a length of 1?
    pts=len(rows)
    
    t_T_P_AMol = [rows[i][0:4] for i in range(0,pts)] 
    t_T_P_AMol = reduce_resolution(t_T_P_AMol,100) 
    
    t_SMol = [[rows[i][0],rows[i][4]] for i in range(0,pts)]  
    t_SMol = reduce_resolution(t_SMol,100) 
    
    t_AllSpecieSens = [[rows[i][0],rows[i][5]] for i in range(0,pts)]
    t_AllSpecieSens = reduce_resolution(t_AllSpecieSens,100)
    
    All_time_Sens.append([x for x in t_AllSpecieSens])
    All_tTP_AMo.append([x for x in t_T_P_AMol])
    return t_T_P_AMol, t_SMol, t_AllSpecieSens


def reduce_resolution(mylist,maxlength):
    """ Reduces the number of elements in a list if the list length is greater
    than the specificed maximum. The function saves every nth term from the list 
    to a new list, where n is calculated to be an optimal sequential number to 
    evenly represent the original list. The resulting number of elements will be 
    less than or equal to the speicified maximum. 
    
    Parameters 
    ----------
    mylist : list
        List for which reduced resolution is desired.
    maxlength : int
        Integer greater than zero which sets the maximum number of elements 
        desired in the list. 
    
    Returns
    ---------
    reduced : list 
        List with a reduced resolution less than or equal to the maxlength. 
    
    """

def reduce_resolution(rows1,max_pts):
   
    reduced =[]
    length=len(mylist)
    if length > maxlength:
        nth = math.ceil(length/maxlength)
        reduced = mylist[::nth]
    else:

        reduced = mylist
    return reduced        

    reduced = rows1
    return reduced      
  

#def reduce_resolution_avg(rowsl,max_pts):
#    
#    reduced=[]
#    pts=len(rows1)
#    if pts > max_pts:
#        nth = math.ceil(pts/max_pts)
#        for i in range(0,pts-nth):
#            for j in range(0,rows1[0])
#        reduced = rows1[::nth]
#    else:
#        reduced = rows1
#    return reduced        

def mole_fractions():
    """Function which zeros mole fractions values.
    
    This function is checking if the mole fractions from the t_SMol
    list are above a predefined ppm (part per million). In this case, 
    if the mole fractions are not above one ppm then the values get
    replaced with zeros. This will facilitate future analysis and 
    ranking of the sensitivities to reactions of interest.
    
    Parameters
    ----------
    SpecificSpecies : list
        List created by user to identify Species of interest.
    molfrac_time : list
        List which contains all of the time steps for the simulation
        being considered.
    molfrac : list
        List which contains all of the mole fractions relating to each
        time step in the 'molfrac_time' list for 'SpecificSpecies.'
    ppm : float
        Floating point number which will define how large the mole
        fractions should be, in ppm.
        
    Dummy Variables
    ---------------
    row : list
        Dummy variable used to replace actual mole fraction values with 
        zero if they are not above the predefined 'ppm' value
        
    Appends
    -------
    molfrac_conditions : list
        List which contains current simulation mixture, Temperature,
        Pressure, mole fractions, and time steps. This list contains
        the updated mole fraction values after checking if they are
        above the predefined ppm value. If they are above the ppm value,
        they are left alone, if they are not they are replaced with 0.
    MoleFraction : list
        This is a list of lists. Each list is identical to the
        information in 'molfrac_conditions' and each list corresponds   
        to a different simulation or initial conditions.
    
    """
    molefrac_time=np.array([x[0] for x in t_SMol])
    molfrac = np.absolute(np.array([x[1] for x in t_SMol]))  #specific specie moles frac
    ppm= 1/1000000  #one ppm
    
    molfrac_conditions = [None]*(len(SpecificSpecies)+2)
    molfrac_conditions[0:2] = [mix,temp, pressure]
    
    for i in range(0,len(SpecificSpecies)):
        row=np.zeros((len(molfrac),2))
        for j in range(0,len(molfrac)):
            if molfrac[j,i] >= ppm:
                row[j,0] = molfrac[j,i]
                row[j,1] = molefrac_time[j]
            else:
                row[j,0] = 0
                row[j,1] = molefrac_time[j]
        molfrac_conditions[i+3] = row  
    MoleFraction.append([x for x in molfrac_conditions])
    return molfrac_conditions


def specific_sens():
    """ Function which zeroes sensitivity values.
    
    This function iterates through all of the 'SpecificSpecie' 
    'molefrac_conditions', if the mole fractions are above the predefined
    ppm value within the 'mole_fractions' function then the sensitivities
    get saved, if not the sensitivities are replaced with value of 0.
    
    Parameters
    ----------
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural 
        gas combustion, including NO formation and reburn chemistry.
    SpecificSpecies : list
        List created by user to identify Species of interest.
    t_SMol : list
        List of lists containing time as the first element, and the 
        corresponding mole fractions, of the chosen SpecificSpecies, as a
        list for the second element.
    SpecificSpecieNumbers : list
        List containing the index number of species in 'SpecificSpecies'
        list with respect to the used mechanism file. (GRI - Mech 3.0).
    molfrac_conditions : list
        Please refer to 'mole_fractions' function.
    sensitivities : list
        List containing sensitivity values for all species in GRI - Mech
        mechanism file with respect to each reaction within the same 
        mechanism file.
        
    Dummy Variables
    ---------------
    row : list
        Dummy variable used for iteration purposes and for appending
        SpecificSpecieSens.
    molfractions : list
        List which grabs the mole fractions and time steps from the
        molfrac_conditions list.
    MolFractions : list
        List which grabs the mole fraction values from molfractions list
        defined above.
    MolFractionsTime : list
        List which grabs the time step values from the molfractions list
        defined above.
        
    Appends
    -------
    SpecificSpecieSens: list
        List which contains the sensitivities of 'SpecificSpecies.' This
        list is dependent on the 'mole_fractions' function, if the mole
        fractions are not above the predefined ppm value then the 
        sensitivites get replaced with zero, otherwise they are 
        not changed. This list is unique and depends on initial
        conditions or mixture of current simulation. This is a list of
        lists, each list corresponds to a specific time step and each 
        list contains the sensitivities of the 'SpecificSpecies' per each 
        reaction. The size of each list should be 
        len(SpecificSpecies)*len(rxns)
        
    """
    row = [None]*len(SpecificSpecies)*len(rxns)
    molefrac_time=np.array([x[0] for x in t_SMol])
    
    for i in range(0,len(molefrac_time)):
        for k in range(0,len(rxns)):
            for j in range(0,len(SpecificSpecieNumbers)):
                molfractions=molfrac_conditions[3+j]
                MolFractions = molfractions[:,0]
                MolFractionsTime = molfractions[:,1]
                if MolFractions[i] == 0:
                    row[j+k*len(SpecificSpecies)] = 0
                elif MolFractions[i] != 0:
                    row[j+k*len(SpecificSpecies)] = sensitivities[i,SpecificSpecieNumbers[j],k]
        SpecificSpecieSens.append([x for x in row])
    return SpecificSpecieSens


def sensitivity_score():
    """Ratio of sum of all maximum reaction values to number of reactions.

    This function obtains the maximum sensitivity from a complete simulation
    per specie, per reaction, and takes the average of all maximum values
    with respect to the number of reactions within the used mechanism file.

    Parameters
    ----------
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural 
        gas combustion, including NO formation and reburn chemistry.
    SpecificSpecies : list
        List created by user to identify Species of interest.
        
    Dummy Variables
    ---------------
    dataa : ndarray
        N dimensional array which contains a list of lists. Each list 
        corresponds to a different time step within a specific simulation.
        Within each list, there are sensitivity values for each 
        'SpecificSpecies' for each reaction within the used mechanism file.
        The size of each list should be len(SpecificSpecies)*len(rxns)
    row : list
        List which contains Temperature, Pressure and the average of the 
        found maximum sensitivities.
    row1 : list
        List which contains the average time at which all maximum 
        sensitivities were found.
    rxn_maxs : ndarray
        N - dimensional array used to collect maximum sensitivities per
        reaction per species throughout all time steps in current 
        simulation results being analyzed.
    rxn_t : ndarray
        N - dimensional array used to store the time at which the maximum
        sensitivities from 'rxn_maxs' occur in the simulation results.
    
            
    Appends
    -------
    score_T_P_MaxSensavg : list
       A list in which elements 0 and 1 corresponds to Temperature,
       Pressure, and elements 2-6 correspond to the ratio of the sum
       of all maximum sensitivites to the number of reactions in the 
       mechanism file being used.
    scoretimes : list
        A list which relates time to elements 2-len(SpecificSpecies) of
        score_T_P_MaxSensavg.
    
    """
    dataa=np.absolute(np.array([x[:] for x in SpecificSpecieSens])) 
    row = [None]*(len(SpecificSpecies)+2)
    row1 = [None]*(len(SpecificSpecies))
    
    row[0] = temp                 #in K
    row[1] = pressure*101.325     #in kPa 
    
    rxn_maxs = np.zeros(len(rxns))
    rxn_t = np.zeros(len(rxns))
    
    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(rxns)):
            rxn_maxs[j] = max(dataa[:,i+j*len(SpecificSpecies)])
            rxn_t[j]=senstime[np.argmax(dataa[:,i+j*len(SpecificSpecies)])]
        row[i+2] = sum(rxn_maxs)/len(rxns) #avg of max sensitivity per reactions
        row1[i] = sum(rxn_t)/len(rxns)
    score_T_P_MaxSensavg.append([x for x in row]) #T,P,avg between maximum sens per reaction
    scoretimes.append([x for x in row1]) #what does this populate?
   
    
def sensitivity_score2(): # finds the maximum sensitivity per reaction per specie
    dataa2=np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    ss = [None]*(len(SpecificSpecies)+2)  #ss-sensititvity score 
    ss_time = [None]*(len(SpecificSpecies))
    max_rxn = [None]*(len(SpecificSpecies))
    
    ss[0] = temp
    ss[1] = pressure*101.325     #in kPa for the plot
    
    max_sens = np.zeros(len(rxns))
    rxn_t = np.zeros(len(rxns))
    
    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(rxns)):
            max_sens[j] = max(dataa2[:,i+j*len(SpecificSpecies)])
            rxn_t[j]=senstime[np.argmax(dataa2[:,i+j*len(SpecificSpecies)])]
        ss[i+2] = max(max_sens)                     #maximum sensitivity
        ss_time[i] = rxn_t[np.argmax(max_sens)]     #time of max sensitivity
        max_rxn[i] =rxns[np.argmax(max_sens)]       #rxn with max sensitivity
    score2_T_P_MaxSens.append([x for x in ss])
    score2times.append([x for x in ss_time])
    score2_Max_sens_rxn.append([x for x in max_rxn])    
    
    
def rank_all(SpecificSpecieSens):
    """ Creates a list of reaction rankings for all timesteps and all species.
        Needs work. Understanding output rank_mat: The lists cycle though the 
        species list, and after every len(species) its a new time_step."""
    data3 = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    step = [None]*len(rxns)
    for i in range(0,len(senstime)):
        for j in range(0,len(SpecificSpecies)):
            for k in range(0,len(rxns)):
                step[k] = data3[i,j+k*len(SpecificSpecies)] #creates list of specie sensitivities per rxns, per time, per specie
            ranking = ranking_per_step(step)
            all_ranks.append([x for x in ranking])
    All_Ranks.append([x for x in all_ranks])
    return all_ranks

def ranking_per_step(sval):
    """ Assigns a rank to each reaction sensitivity per time-step and species""" 
    sval = np.array(sval) # sval = step from above per time, per specie
    order = np.argsort(sval) # orders sval from smallest to largeest by listing indices which sval elements should be arranged in
    ranks = ranks_possible(sval)
    step_rank = [None]*len(sval)
    for j in range (0,len(sval)):
        step_rank[order[j]] = ranks[j]
    return step_rank
    
def ranks_possible(sval):
    """ Determines how many ranks there will be for a given time-step, since 
    sensitivities can tie if they are equal"""
    sval = np.sort(sval) #orders elements from smallest to largest
    sval = np.flip(sval,axis=0) #re orders from largest to smallest
    ranks = [None]*len(sval)
    if sval[0] == 0:
        ranks[0] = len(sval)
    else:
        ranks[0] = 1
    for i in range (1,len(sval)):
        if sval[i] == 0:  
            ranks[i] = len(sval)
        else:
            if sval[i] == sval[i-1]:
                ranks[i] = ranks[i-1]
            else:
                ranks[i] = ranks[i-1]+1
    ranks = ranks[::-1]             #Flips the list 
    return ranks


def interested_rxn_ranking(all_ranks):
    """ Takes the rankings of all reactions and returns a list of the ranking 
    of the reactions of interest for each species for each condition and the 
    initial inputs for that condition""" 
    all_ranks = np.array(all_ranks)
    indices = rxns_interested_indices(rxns,interested_rxns) 
    time_steps = len(senstime)
    species_sens = [None]*time_steps
    rank_across_spec =[None]*len(indices)
    ranks_per_cond = []
    rank = [None]*5
    rank[0] = phi
    rank[1] = methane
    rank[2] = temp
    rank[3] = pressure*101.325
    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(indices)):
            sense = all_ranks[:,indices[j]] 
            for k in range(0,time_steps):
                species_sens[k] = sense[i+k*len(SpecificSpecies)]
            rank_across_spec[j] = min(species_sens)
        ranks_per_cond.append([x for x in rank_across_spec])
        rank_plot.append([x for x in rank_across_spec])
    rank[4] = np.array(ranks_per_cond)
    rank_score.append([x for x in rank])
    return rank_score, rank_plot

def rxns_interested_indices(rxns,interested_rxns):
    """ Creates an array of indicies of the reactions of interest with respect
    to the list of all reactions"""
    rxns = rxns.tolist()
    rxn_indx = []
    for rxn in interested_rxns:
        ind = rxns.index(rxn)
        rxn_indx.append(ind)
    return rxn_indx
   
    
def ranking_plots(rank_score,rank_plot):
    
    Phi_ary=[None]*len(rank_score)
    CH4_ary=[None]*len(rank_score)
    T_ary=[None]*len(rank_score)
    P_ary=[None]*len(rank_score)
    for cond in range(0,len(rank_score)):    
        Phi_ary[cond]=rank_score[cond][0]
        CH4_ary[cond]=rank_score[cond][1]
        T_ary[cond]=rank_score[cond][2]
        P_ary[cond]=rank_score[cond][3]    
    Phi_ary=np.array(Phi_ary)
    CH4_ary=np.array(CH4_ary)
    T_ary=np.array(T_ary)
    P_ary=np.array(P_ary)
    condition_num = int(len(rank_plot)/len(SpecificSpecies))
    for rxn in interested_rxns:
        ranks = [None]*condition_num
        fig1 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax1 = fig1.add_subplot(111)
        ax1.grid()
        ax1.set_ylabel('Ranking')
        ax1.set_xlabel('Pressure [kPa]')
        ax1.set_title('Rxn '+str(rxn)+' Ranking vs Presure')
        fig2 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax2 = fig2.add_subplot(111)
        ax2.grid()
        ax2.set_ylabel('Ranking')
        ax2.set_xlabel('Temperature [K]')
        ax2.set_title('Rxn '+str(rxn)+' Ranking vs Temperature')
        fig3 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax3 = fig3.add_subplot(111)
        ax3.grid()
        ax3.set_ylabel('Ranking')
        ax3.set_xlabel('Methane.X')
        ax3.set_title('Rxn '+str(rxn)+' Ranking vs Methane Mole Fraction')
        fig4 = plt.figure()
        plt.ylim(len(rxns)+1,0) 
        ax4 = fig4.add_subplot(111)
        ax4.grid()
        ax4.set_ylabel('Ranking')
        ax4.set_xlabel('Equivalence Ratio')
        ax4.set_title('Rxn '+str(rxn)+' Ranking vs Equivalence Ratio')
        marker=['o','s','v','^','>']
        color=['b','darkorange','g','r','blueviolet']
        for species, mrk, c in zip(SpecificSpecies,marker,color):
            for i in range(0,condition_num):
                ranks[i] = rank_plot[i*len(SpecificSpecies)+
                                       SpecificSpecies.index(species)][interested_rxns.index(rxn)]
            ranks = np.array(ranks)
            ax1.scatter(P_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig1.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig1.tight_layout()
            fig1.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Pressure.png',bbox_inches = 'tight')
            ax2.scatter(T_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig2.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig2.tight_layout()
            fig2.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Tempurature.png',bbox_inches = 'tight')
            ax3.scatter(CH4_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig3.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig3.tight_layout()
            fig3.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Methane%.png',bbox_inches = 'tight')
            ax4.scatter(Phi_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig4.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig4.tight_layout()
            fig4.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Phi.png',bbox_inches = 'tight')


def plot_directory():
    """ Creates directories inside the working directory for all the plots in 
        the script. current_directory/Plots/date/file_number/plot_type/file """ 
        
    file_num = 1
    current_dir = os.getcwd()
    ranks_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Ranking')
    while os.path.exists(ranks_dir) == True: 
        file_num = file_num+1
        ranks_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Ranking')    
    if not os.path.exists(ranks_dir):
        os.makedirs(ranks_dir)
    score2_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Score2')
    if not os.path.exists(score2_dir):
        os.makedirs(score2_dir)   
    score_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Score')
    if not os.path.exists(score_dir):
        os.makedirs(score_dir) 
    data_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)        
    return file_num

def progress_track(CH4,methane,Phi,phi):
    """Prints the progress status of the code runtime as a percentage based on the 
    iteration of the CH4 parameter."""
    CH4_list= CH4.tolist()
    Phi_list= Phi.tolist()
    CH4_iteration= CH4_list.index(methane)
    Phi_iteration=Phi_list.index(phi)
    progress= str(round((CH4_iteration+len(CH4)*Phi_iteration)/(len(Phi)*len(CH4))*100))
    print(progress+'%')

        
if __name__ == "__main__":
    #first scoring criteria
    score_T_P_MaxSensavg = []    #T, P, avg of max values between all rxns per specific specie 
    scoretimes = []
    rank_score=[]
    rank_plot=[]
    All_Ranks=[]
    #second scoring criteria
    score2_T_P_MaxSens=[]
    score2times=[]
    score2_Max_sens_rxn=[]  
    #mole/massfraction function
    MoleFraction = []
    All_time_Sens=[]
    All_tTP_AMo=[]
    dt = datetime.datetime.now()
    date = str(dt.month)+'-'+str(dt.day)
    tic=time.time()
    #create gas object
    gas1 = cantera.Solution('grimech30.cti') #53 species, 325 reactions
    names=gas1.species_names
    SpecificSpecieNumbers=[names.index(speci) for speci in SpecificSpecies]
    for phi in Phi:
        for methane in CH4:
            oxygen = 2/phi*methane
            n2 = 1 - oxygen - methane
            if n2<0:
               continue 
            mix={'CH4':methane, 'O2':oxygen, 'N2': n2}
            progress_track(CH4,methane,Phi,phi)
            
            for temp in T:
                for pressure in P:
                    #set state 1atm <= P <= 30atm, 1atm = 101325 Pa
                    gas1.TP = temp, pressure*101325           
                    #set gas compsotion; X - mole fraction;  Y - mass fraction
                    #gas1.X= 'CH4:0.05, O2:0.2, N2:0.75'
                    gas1.X= mix
                    #create 0D - simulation objects
    #               reac = cantera.IdealGasReactor(gas1)
    #               reac = cantera.ConstPressureReactor(gas1)
                    reac = cantera.IdealGasConstPressureReactor(gas1)
                    sim = cantera.ReactorNet([reac])
                    #initialize parameters
                    SpecificSpecieSens = []   #list of sens of interest per rxn
                    all_ranks = []
                    [t_T_P_AMol, t_SMol, t_AllSpecieSens]= reac_sens()
                    #counter=counter+1
                    senstime = np.array([x[0] for x in t_AllSpecieSens])   #times for all specie sensitivities
                    sensitivities = np.absolute(np.array([x[1] for x in t_AllSpecieSens])) #sensitivities
                    molfrac_conditions=mole_fractions()
                    specific_sens()
                    sensitivity_score()
                    sensitivity_score2()
                    rank_all(SpecificSpecieSens)
                    interested_rxn_ranking(all_ranks) 
    toc=time.time()
    duration = toc-tic
    numbr_mix = len(Phi)
    tp_combs = len(T)*len(P)
    print(numbr_mix, 'mixture with',tp_combs, 'Temperature/Pressure combinations took', duration, 'seconds.')
#    file_num = plot_directory()
#    ranking_plots(rank_score,rank_plot)

#File saving 
#f1=MoleFraction
#with open('Plots/'+date+'/'+str(file_num)+'/Data/MoleFraction.csv', "wb") as fp1:   #Pickling
#    pickle.dump(f1,fp1) 
#fp1.close()    
#f2=SpecificSpecies
#with open('Plots/'+date+'/'+str(file_num)+'/Data/SpecificSpecies.csv', "wb") as fp2:   #Pickling
#    pickle.dump(f2,fp2) 
#fp2.close() 
#f3=All_tTP_SMo_AMo_SMa_AMa
#with open('Plots/'+date+'/'+str(file_num)+'/Data/All_tTP_SMo_AMo_SMa_AMa.csv', "wb") as fp3:   #Pickling
#    pickle.dump(f3,fp3) 
#fp3.close()
#f4=All_time_Sens
#with open('Plots/'+date+'/'+str(file_num)+'/Data/All_time_Sens.csv', "wb") as fp4:   #Pickling
#    pickle.dump(f4,fp4) 
#fp4.close()     
#f5=T
#with open('Plots/'+date+'/'+str(file_num)+'/Data/T.csv', "wb") as fp5:   #Pickling
#    pickle.dump(f5,fp5) 
#fp5.close()
#f6=P
#with open('Plots/'+date+'/'+str(file_num)+'/Data/P.csv', "wb") as fp6:   #Pickling
#    pickle.dump(f6,fp6) 
#fp6.close() 
#f7=Phi
#with open('Plots/'+date+'/'+str(file_num)+'/Data/Phi.csv', "wb") as fp7:   #Pickling
#    pickle.dump(f7,fp7) 
#fp7.close()
#f8=CH4
#with open('Plots/'+date+'/'+str(file_num)+'/Data/CH4.csv', "wb") as fp8:   #Pickling
#    pickle.dump(f8,fp8) 
#fp8.close()   
#f9=All_Ranks
#with open('Plots/'+date+'/'+str(file_num)+'/Data/All_ranks.csv', "wb") as fp9:   #Pickling
#    pickle.dump(f9,fp9) 
#fp9.close()
#f10=rank_score
#with open('Plots/'+date+'/'+str(file_num)+'/Data/rank_score.csv', "wb") as fp10:   #Pickling
#    pickle.dump(f10,fp10) 
#fp10.close()
#f11=rank_plot
#with open('Plots/'+date+'/'+str(file_num)+'/Data/rank_plot.csv', "wb") as fp11:   #Pickling
#    pickle.dump(f11,fp11) 
#fp11.close()

#plotting

#'%matplotlib qt' --> type in console to seperate plots from console
#sens_score2 plots
#score2_T_P_MaxSens = np.array([x[:] for x in score2_T_P_MaxSens])   #scores: avg of max sens per rxn, per specie
#score2times = np.array([x[:] for x in score2times])   # times: time avg of max sens per rxn per specie
#score2_Max_sens_rxn = np.array([x[:] for x in score2_Max_sens_rxn]) 
#plt.figure('H20')
#plt.subplot(211)
#plt.scatter(score2_Max_sens_rxn[:,0],score2_T_P_MaxSens[:,2],label='H2O')
#plt.ylabel('Sensitivity')
#plt.xlabel('Reaction')
#plt.title('Reaction/Time vs. Sensitivity for H2O')
#plt.grid(axis='both')
#plt.subplot(212)
#plt.scatter(score2times[:,0],score2_T_P_MaxSens[:,2],label='H2O')
#plt.ylabel('Sensitivity')
#plt.xlabel('Time')
##plt.legend()
#plt.minorticks_on()
#plt.grid(axis='both') 
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score2/H20.png')
#
#plt.figure('CH4')
#plt.subplot(211)
#plt.scatter(score2_Max_sens_rxn[:,1],score2_T_P_MaxSens[:,3],label='CH4')
#plt.ylabel('Sensitivity')
#plt.xlabel('Reaction')
#plt.title('Reaction/Time vs. Sensitivity for CH4')
#plt.grid(axis='both')
#plt.subplot(212)
#plt.scatter(score2times[:,1],score2_T_P_MaxSens[:,3],label='CH4')
#plt.ylabel('Sensitivity')
#plt.xlabel('Time')
##plt.legend()
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score2/CH4.png')
#
#plt.figure('CO')
#plt.subplot(211)
#plt.scatter(score2_Max_sens_rxn[:,2],score2_T_P_MaxSens[:,4],label='CO')
#plt.ylabel('Sensitivity')
#plt.xlabel('Reaction')
#plt.title('Reaction/Time vs. Sensitivity for CO')
#plt.grid(axis='both')
#plt.subplot(212)
#plt.scatter(score2times[:,2],score2_T_P_MaxSens[:,4],label='CO')
#plt.ylabel('Sensitivity')
#plt.xlabel('Time')
##plt.legend()
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score2/C0.png') 
#
#plt.figure('CO2')
#plt.subplot(211)
#plt.scatter(score2_Max_sens_rxn[:,3],score2_T_P_MaxSens[:,5],label='CO2')
#plt.ylabel('Sensitivity')
#plt.xlabel('Reaction')
#plt.title('Reaction/Time vs. Sensitivity for CO2')
#plt.grid(axis='both')
#plt.subplot(212)
#plt.scatter(score2times[:,3],score2_T_P_MaxSens[:,5],label='CO2')
#plt.ylabel('Sensitivity')
#plt.xlabel('Time')
##plt.legend()
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score2/CO2.png')
#
#plt.figure('NO')
#plt.subplot(211)
#plt.scatter(score2_Max_sens_rxn[:,4],score2_T_P_MaxSens[:,6],label='NO')
#plt.ylabel('Sensitivity')
#plt.xlabel('Reaction')
#plt.title('Reaction/Time vs. Sensitivity for NO')
#plt.grid(axis='both')
#plt.subplot(212)
#plt.scatter(score2times[:,4],score2_T_P_MaxSens[:,6],label='NO')
#plt.ylabel('Sensitivity')
#plt.xlabel('Time')
##plt.legend()
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score2/NO.png')
#
#
##sens_score plots
#score_T_P_MaxSensavg = np.array([x[:] for x in score_T_P_MaxSensavg])   #scores: avg of max sens per rxn, per specie
#scoretimes = np.array([x[:] for x in scoretimes])   # times: time avg of max sens per rxn per specie
#plt.figure('Score v. T')
#plt.scatter(score_T_P_MaxSensavg[:,0],score_T_P_MaxSensavg[:,2],label='H2O')
#plt.scatter(score_T_P_MaxSensavg[:,0],score_T_P_MaxSensavg[:,3],label='CH4')
#plt.scatter(score_T_P_MaxSensavg[:,0],score_T_P_MaxSensavg[:,4],label='CO')
#plt.scatter(score_T_P_MaxSensavg[:,0],score_T_P_MaxSensavg[:,5],label='CO2')
#plt.scatter(score_T_P_MaxSensavg[:,0],score_T_P_MaxSensavg[:,6],label='NO')
#plt.ylabel('Avg of all max sensitivities')
#plt.xlabel('Temperature [K]')
#plt.title('Score vs Temperature')
#plt.legend(bbox_to_anchor=(1.02, 1), loc=2)
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score/Temperature.png')
#
#plt.figure('Score v. P')
#plt.scatter(score_T_P_MaxSensavg[:,1],score_T_P_MaxSensavg[:,2],label='H2O')
#plt.scatter(score_T_P_MaxSensavg[:,1],score_T_P_MaxSensavg[:,3],label='CH4')
#plt.scatter(score_T_P_MaxSensavg[:,1],score_T_P_MaxSensavg[:,4],label='CO')
#plt.scatter(score_T_P_MaxSensavg[:,1],score_T_P_MaxSensavg[:,5],label='CO2')
#plt.scatter(score_T_P_MaxSensavg[:,1],score_T_P_MaxSensavg[:,6],label='NO')
#plt.ylabel('Avg of all max sensitivities')
#plt.xlabel('Pressure [kPa]')
#plt.title('Score vs Pressure')
#plt.legend(bbox_to_anchor=(1.02, 1), loc=2)
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score/Pressure.png')
#
#plt.figure('Score v. time')
#plt.scatter(scoretimes[:,0],score_T_P_MaxSensavg[:,2],label='H2O')
#plt.scatter(scoretimes[:,1],score_T_P_MaxSensavg[:,3],label='CH4')
#plt.scatter(scoretimes[:,2],score_T_P_MaxSensavg[:,4],label='CO')
#plt.scatter(scoretimes[:,3],score_T_P_MaxSensavg[:,5],label='CO2')
#plt.scatter(scoretimes[:,4],score_T_P_MaxSensavg[:,6],label='NO')
#plt.ylabel('Avg of all max sensitivities')
#plt.xlabel('Time [s]')
#plt.title('Score vs Time')
#plt.legend(bbox_to_anchor=(1.02, 1), loc=2)
#plt.minorticks_on()
#plt.grid(axis='both')
#plt.tight_layout()
#plt.savefig('Plots/'+date+'/'+str(file_num)+'/Score/Time.png') 