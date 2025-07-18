#!/usr/bin/env python

## \file mmg.py
#  \brief python script for running mesh adaptation using the MMG Inria library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os, shutil, copy, time

from SU2 import io as su2io
from SU2.adap.tools import *
from SU2.adap.interface import run_command, call_mmg
from SU2.run.interface import CFD as SU2_CFD
from CDriver import CDriver

def mmg(config):
    """
    Runs the a mesh adaptation loop with the MMG library.

    Inputs:
        config - an SU2 config object
    """

    #--- Check config options related to mesh adaptation

    pyadap_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 
                      'ADAP_FLOW_ITER', 'ADAP_ADJ_ITER', 'ADAP_CFL', 'ADAP_HAUSD' ]
    required_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_HMAX', 'ADAP_HMIN', 'MESH_FILENAME', 
                        'RESTART_SOL', 'MESH_OUT_FILENAME' ]

    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise AttributeError(err)
    
    #--- NEMO solver check
    if 'NEMO' in config.SOLVER:
        nemo = True
    else:
        nemo = False

    #--- Print adap options

    print(print_adap_options(config))

    #--- Target mesh sizes and subiterations at each size

    mesh_sizes = get_mesh_sizes(config)
    sub_iter   = get_sub_iterations(config)

    if len(mesh_sizes) != len(sub_iter):
        raise ValueError(f'Inconsistent number of mesh sizes and sub-iterations. {len(mesh_sizes)} mesh sizes and {len(sub_iter)} sub-iterations provided.')

    #--- Solver iterations/ residual reduction param for each size level

    flow_iter = get_flow_iter(config)
    flow_cfl  = get_flow_cfl(config)

    adap_sensors = get_adap_sensors(config)
    sensor_avail = ['MACH', 'PRESSURE', 'TEMPERATURE', 'ENERGY', 'DENSITY', 'TOTALPRESSURE']

    for sensor in adap_sensors:
        if sensor not in sensor_avail:
            raise ValueError(f'Unknown adaptation sensor {sensor}. Available options are {sensor_avail}.')
        
    sensor = 'Pressure'

    #--- Change current directory

    warn = True
    base_dir = os.getcwd()
    adap_dir = './adap'
    flow_dir = './Flows'
    surfflow_dir = './Flows_surf'

    if os.path.exists(adap_dir):
        print('./adap exists. Removing old mesh adaptation in 10s.')
        if warn : time.sleep(10)
        shutil.rmtree(adap_dir)
        print(f'The {adap_dir} folder was deleted.')
    
    if os.path.exists(flow_dir):
        print(flow_dir+' exists. Removing old Flows in 10s.')
        shutil.rmtree(flow_dir)
        print(f'The {flow_dir} folder was deleted.')
    
    if os.path.exists(surfflow_dir):
        print(surfflow_dir+' exists. Removing old Surface Flows in 10s.')
        shutil.rmtree(surfflow_dir)
        print(f'The {surfflow_dir} folder was deleted.')


    os.makedirs(flow_dir)
    os.makedirs(surfflow_dir)
    run_command(f'ln -s ../adap/ite{0}/'+config.VOLUME_FILENAME+'.vtu ./Flows/'+config.VOLUME_FILENAME+'_'+str(0).zfill(5)+'.vtu')
    run_command(f'ln -s ../adap/ite{0}/'+config.SURFACE_FILENAME+'.vtu ./Flows_surf/'+config.SURFACE_FILENAME+'_'+str(0).zfill(5)+'.vtu')
    dir = f'{adap_dir}/ite0'
    os.makedirs(dir)
    os.chdir(dir)
    os.symlink(os.path.join(base_dir, config.MESH_FILENAME), config.MESH_FILENAME)

    meshfil = config['MESH_FILENAME']

    #--- Format of history file

    history_format = config.TABULAR_FORMAT
    if (history_format == 'TECPLOT'):
        history_filename = os.path.join(base_dir, 'history_adap.dat')
    else:
        history_filename = os.path.join(base_dir, 'history_adap.csv')

    #--- Get mesh dimension

    dim = get_su2_dim(meshfil)
    if ( dim != 2 and dim != 3 ):
        raise ValueError('Wrong dimension number.')

    #--- MMG parameters

    config_mmg = get_mmg_config(config, dim)
    config_mmg['toll'] = 1
    config_mmg['card'] = 800

    #--- Compute initial solution if needed, else link current files

    config_cfd = copy.deepcopy(config)
    for opt in pyadap_options:
        config_cfd.pop(opt, None)

    #--- Check config for filenames if restarting
    restart = config['RESTART_SOL'] == 'YES'
    if restart:
        required_options = ['SOLUTION_FILENAME']
        if not all (opt in config for opt in required_options):
            err = 'RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise ValueError(err)

        os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), config.SOLUTION_FILENAME)

        print('\nInitial CFD solution is provided.')

    else:
        print('\nRunning initial CFD solution.')

    #--- Only allow ASCII restarts for file conversion AGGIUNTA LETTURA BINARIO
    if '.csv' in config.RESTART_FILENAME:
        config_cfd.READ_BINARY_RESTART = 'NO'
        sol_ext_cfd = '.csv'
        config_cfd.OUTPUT_FILES = ['RESTART_CSV','PARAVIEW','SURFACE_PARAVIEW']
    else:
        sol_ext_cfd = '.dat'
        config_cfd.OUTPUT_FILES = ['RESTART','PARAVIEW','SURFACE_PARAVIEW']


    solfil  = f'restart_flow{sol_ext_cfd}'
    set_flow_config_ini(config_cfd, solfil, adap_sensors, mesh_sizes[0])

    try: # run with redirected outputs
        #--- Run a single iteration of the flow if restarting to get history info
        if restart:
            config_cfd.ITER = 1
            config_cfd.RESTART_CFL = 'NO'

        with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

        #--- Set RESTART_SOL=YES for runs after adaptation
        if not nemo and sol_ext_cfd != '.csv':     # SU2 does not interpolate ASCII files for restarts
            config_cfd.RESTART_SOL = 'YES' 
            config_cfd.RESTART_CFL = 'NO'

    except:
        raise

    #--- Check existence of initial mesh, solution

    required_files = [meshfil, solfil]

    if not all (os.path.exists(fil) for fil in required_files):
        err = "Can't find the following files:\n"
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise Exception(err)

    #--- Start adaptive loop

    global_iter = 0

    #--- Print convergence history

    npoin = get_su2_npoin(meshfil)
    plot_results(history_format, history_filename, global_iter, npoin)

    print('\nStarting mesh adaptation process.\n')

    nSiz = len(mesh_sizes)
    for iSiz in range(nSiz):
        nSub = int(sub_iter[iSiz])
        for iSub in range(nSub):
            
            global_iter += 1

            os.symlink(f'../adap/ite{global_iter}/'+config_cfd.VOLUME_FILENAME+'.vtu', '../../Flows/'+config_cfd.VOLUME_FILENAME+'_'+str(global_iter).zfill(5)+'.vtu')
            os.symlink(f'../adap/ite{global_iter}/'+config_cfd.SURFACE_FILENAME+'.vtu', '../../Flows_surf/'+config_cfd.SURFACE_FILENAME+'_'+str(global_iter).zfill(5)+'.vtu')

            mesh_size = int(mesh_sizes[iSiz])
            if iSub == nSub-1 and iSiz != nSiz-1: mesh_size = int(mesh_sizes[iSiz+1])
            config_mmg['size'] = mesh_size

            #--- Instantiating the driver
            driver = CDriver(sensor, meshfil, solfil,params=config_mmg)

            #--- Initializing the driver (SU2 mesh and sol reading)
            driver.ReadSU2()

            #--- Computing the metric 
            driver.ComputeMetric()

            #--- Writing the MMG mesh, sol and param file
            driver.WriteMedit(meshfil.replace('.su2', '.mesh'),
                              solfil.replace(sol_ext_cfd, '.sol'))


            #--- Adapt mesh with MMG
            meshin = config_cfd['MESH_FILENAME'].replace('.su2','.mesh')
            meshout = config_cfd['MESH_OUT_FILENAME'].replace('.su2','.mesh')

            solfile = config_cfd['RESTART_FILENAME'].replace(sol_ext_cfd,'.sol')
            call_mmg(meshin, meshout, solfile, config_mmg)

            #--- Reading the adapted mesh in MMG format
            driver.ReadMedit(meshout)

            #--- Writing the adapted mesh in SU2 format
            driver.WriteSU2(meshout.replace('.mesh','.su2'))

            #--- Print mesh sizes
            print_adap_table(iSiz, mesh_sizes, iSub, nSub, driver.mesh.meshDict)

            dir = f'./ite{global_iter}'
            os.makedirs(os.path.join('..',dir))
            os.chdir(os.path.join('..',dir))

            meshfil = config_cfd['MESH_FILENAME']
            solfil_ini  = f'restart_flow{sol_ext_cfd}'

            driver.WriteSU2(meshfil)

            #--- Run su2

            # link to restart if restart file is binary
            if sol_ext_cfd != '.csv':
                os.symlink(os.path.join(f'../ite{global_iter-1}', config.RESTART_FILENAME), solfil_ini)

            try: # run with redirected outputs

                update_flow_config(config_cfd, meshfil, solfil, solfil_ini,
                                        flow_iter[iSiz], flow_cfl[iSiz], adap_sensors, mesh_size)

                with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

                if not os.path.exists(solfil) :
                    raise RuntimeError('SU2_CFD failed.\n')

                #--- Print convergence history

                npoin = get_su2_npoin(meshfil)
                plot_results(history_format, history_filename, global_iter, npoin)

            except:
                raise

            del driver

    #--- Write final files

    # fileconverter = MeshSolConverter()
    # fileconverter.SU2ToMeditMesh(meshfil, meshfil.replace('.su2', '.mesh'))

    # os.rename(solfil, os.path.join(base_dir, config.RESTART_FILENAME))
    # os.rename(meshfil, os.path.join(base_dir, config.MESH_OUT_FILENAME))

    pad_nul = ' '*15
    print('\nMesh adaptation successfully ended.')
    print(f'Results files: {config.MESH_OUT_FILENAME}\n{pad_nul}{config.RESTART_FILENAME}')
