##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf.py                                                          #
#    Nexus interface to the PWSCF simulation code.                   #
#                                                                    #
#  Content summary:                                                  #
#    Pwscf                                                           #
#      Simulation class for PWSCF.                                   #
#                                                                    #
#    generate_pwscf                                                  #
#      User-facing function to create Pwscf simulation objects.      #
#                                                                    #
#====================================================================#


import os
from numpy import array
from generic import obj
from physical_system import PhysicalSystem
from simulation import Simulation
from pwscf_input import PwscfInput,generate_pwscf_input
from pwscf_analyzer import PwscfAnalyzer
from execute import execute
from debug import ci


unique_vdw_functionals = [
    'optb86b-vdw',
    'vdw-df3', # optB88+vdW
    'vdw-df',
    'vdw-df2',
    'rev-vdw-df2',
    'vdw-df-c09',
    'vdw-df2-c09',
    'rvv10',
    ]
repeat_vdw_functionals = [
    'vdw-df4', # 'optB86b-vdW'
    ]
unique_functionals = [
    'revpbe',
    'pw86pbe',
    'b86bpbe',
    'pbe0',
    'hse',
    'gaup',
    'pbesol',
    'pbeq2d',
    'optbk88',
    'optb86b',
    'pbe',
    'wc',
    'b3lyp',
    'pbc',
    'bp',
    'pw91',
    'hcth',
    'olyp',
    'tpss',
    'oep',
    'hf',
    'blyp',
    'lda',
    'sogga',
    'm06l',
    'ev93',
    ]+unique_vdw_functionals
repeat_functionals = [
    'q2d', # pbeq2d
    'pz', # lda
    ]+repeat_vdw_functionals

vdw_functionals     = set(unique_vdw_functionals+repeat_vdw_functionals)
allowed_functionals = set(unique_functionals+repeat_functionals)



class Pwscf(Simulation):
    input_type = PwscfInput
    analyzer_type = PwscfAnalyzer
    generic_identifier = 'pwscf'
    application = 'pw.x'
    application_properties = set(['serial','mpi'])
    application_results    = set(['charge_density','orbitals','structure','restart'])

    supports_restarts = True # supports restartable, but not force restart yet

    vdw_table = None

    @staticmethod
    def settings(vdw_table=None):
        # van der Waals family of functional require the vdW table generated by
        #  generate_vdW_kernel_table.x: specify 'vdw_table' in settings
        Pwscf.vdw_table = vdw_table
    #end def settings              

    @staticmethod
    def restore_default_settings():
        Pwscf.vdw_table = None
    #end def restore_default_settings


    #def propagate_identifier(self):
    #    self.input.control.prefix = self.identifier        
    ##end def propagate_identifier


    def __init__(self,**sim_args):
        group_atoms   = sim_args.pop('group_atoms',False)
        sync_from_scf = sim_args.pop('sync_from_scf',True)
        Simulation.__init__(self,**sim_args)
        self.sync_from_scf = False
        calc = self.input.control.get('calculation',None)
        if calc=='nscf':
            self.sync_from_scf = sync_from_scf
        #end if
        if group_atoms and isinstance(self.system,PhysicalSystem):
            self.warn('requested grouping by atomic species, but pwscf does not group atoms anymore!')
            #self.system.structure.group_atoms()
        #end if
    #end def __init__


    def write_prep(self):
        #make sure the output directory exists
        outdir = os.path.join(self.locdir,self.input.control.outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        #end if
        #copy over vdw_table for vdW-DF functional
        if self.path_exists('input/system/input_dft'):
            functional = self.input.system.input_dft.lower()
            if '+' not in functional and functional not in allowed_functionals:
                self.warn('functional "{0}" is unknown to pwscf'.format(functional))
            #end if
            if functional in vdw_functionals:
                if self.vdw_table is None:
                    self.error('attempting to run vdW functional "{0}", but vdw_table is missing\nplease provide path to table file via "vdw_table" parameter in settings'.format(functional))
                #end if
                cd_rel = os.path.relpath(self.vdw_table,self.locdir)
                # copy instead of link to vdw_table to avoid file-lock from multiple pw.x instances
                cp_cmd = 'cd '+self.locdir+';cp '+cd_rel+' .'
                os.system(cp_cmd)
            #end if
        #end if          
    #end def write_prep


    def check_result(self,result_name,sim):
        input = self.input
        control = input.control
        if result_name=='charge_density' or result_name=='restart':
            calculating_result = True
        elif result_name=='orbitals':
            calculating_result = 'calculation' not in control or 'scf' in control.calculation.lower()
        elif result_name=='structure':
            calculating_result = 'calculation' in control and 'relax' in control.calculation.lower()
        else:
            calculating_result = False
        #end if
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()        
        input = self.input
        control = input.control
        prefix = 'pwscf'
        outdir = './'
        if 'prefix' in control:
            prefix = control.prefix
        #end if
        if 'outdir' in control:
            outdir = control.outdir
        #end if
        if outdir.startswith('./'):
            outdir = outdir[2:]
        #end if
        if result_name=='charge_density' or result_name=='restart':
            result.locdir   = self.locdir
            result.outdir   = os.path.join(self.locdir,outdir)
            result.location = os.path.join(self.locdir,outdir,prefix+'.save','charge-density.dat')
            result.spin_location = os.path.join(self.locdir,outdir,prefix+'.save','spin-polarization.dat')
        elif result_name=='orbitals':
            result.location = os.path.join(self.locdir,outdir,prefix+'.wfc1')
        elif result_name=='structure':
            pa = self.load_analyzer_image()
            structs = pa.structures
            struct  = structs[len(structs)-1]
            pos   = struct.positions
            atoms = struct.atoms
            if 'celldm(1)' in self.input.system:
                scale = self.input.system['celldm(1)']
            else:
                scale = 1.0
            #end if
            pos   = scale*array(pos)
            
            structure = self.system.structure.copy()
            structure.change_units('B')
            structure.pos = pos
            structure.set_elem(atoms)
            if 'axes' in struct:
                structure.axes = struct.axes.copy()
            #end if
            result.structure = structure
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        if result_name=='charge_density':
            c = self.input.control
            res_path = os.path.abspath(result.locdir)
            loc_path = os.path.abspath(self.locdir)
            if res_path==loc_path:
                None # don't need to do anything if in same directory
            elif self.sync_from_scf: # rsync output into nscf dir
                outdir = os.path.join(self.locdir,c.outdir)
                command = 'rsync -av {0}/* {1}/'.format(result.outdir,outdir)
                if not os.path.exists(outdir):
                    print('Running rsync for the {} directory. This might take a while.'.format(outdir))
                    os.makedirs(outdir)
                #end if
                sync_record = os.path.join(outdir,'nexus_sync_record')
                if not os.path.exists(sync_record):
                    execute(command)
                    f = open(sync_record,'w')
                    f.write('\n')
                    f.close()
                #end if
            else: # attempt to use symbolic links instead
                link_loc = os.path.join(self.locdir,c.outdir,c.prefix+'.save')
                cd_loc = result.location
                cd_rel = os.path.relpath(cd_loc,link_loc)
                sp_loc = result.spin_location
                sp_rel = os.path.relpath(sp_loc,link_loc)
                cwd = os.getcwd()
                if not os.path.exists(link_loc):
                    os.makedirs(link_loc)
                #end if
                os.chdir(link_loc)
                os.system('ln -s '+cd_rel+' charge-density.dat')
                os.system('ln -s '+sp_rel+' spin-polarization.dat')
                os.chdir(cwd)
            #end if
        elif result_name=='structure':
            relstruct = result.structure.copy()
            relstruct.change_units('B')
            self.system.structure = relstruct
            self.system.remove_folded()

            input = self.input
            preserve_kp = 'k_points' in input and 'specifier' in input.k_points and (input.k_points.specifier=='automatic' or input.k_points.specifier=='gamma')
            if preserve_kp:
                kp = input.k_points.copy()
            #end if
            input.incorporate_system(self.system)
            if preserve_kp:
                input.k_points = kp
            #end if
        elif result_name=='restart':
            c = self.input.control
            if('startingwfc' in self.input.electrons and self.input.electrons.startingwfc != 'file'):
                self.error('Exiting. User has specified startingwfc=\''+self.input.electrons.startingwfc+'\'.\nThis value will be overwritten when incorporating result \'restart\'.\nPlease fix conflict.')
            #end if
            if('startingpot' in self.input.electrons and self.input.electrons.startingpot != 'file'):
                self.error('Exiting. User has specified startingpot=\''+self.input.electrons.startingpot+'\'.\nThis value will be overwritten when incorporating result \'restart\'.\nPlease fix conflict.')
            #end if
            c.restart_mode='restart'
            res_path = os.path.abspath(result.locdir)
            loc_path = os.path.abspath(self.locdir)
            if res_path==loc_path:
                None # don't need to do anything if in same directory
            else: # rsync output into new scf dir
                outdir = os.path.join(self.locdir,c.outdir)
                command = 'rsync -av {0}/* {1}/'.format(result.outdir,outdir)
                if not os.path.exists(outdir):
                    print('Running rsync for the {} directory. This might take a while.'.format(outdir))
                    os.makedirs(outdir)
                #end if
                sync_record = os.path.join(outdir,'nexus_sync_record')
                if not os.path.exists(sync_record):
                    execute(command)
                    f = open(sync_record,'w')
                    f.write('\n')
                    f.close()
                #end if
            #end if
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if        
    #end def incorporate_result


    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        not_converged = 'convergence NOT achieved'  in output
        time_exceeded = 'Maximum CPU time exceeded' in output
        user_stop     = 'Program stopped by user request' in output
        run_finished  = 'JOB DONE' in output
        restartable = not_converged or time_exceeded or user_stop
        restart = run_finished and self.restartable and restartable
        if restart:
            self.save_attempt()
            self.input.control.restart_mode = 'restart'
            self.reset_indicators()
        else:
            error_in_routine = 'Error in routine' in output
            failed = not_converged or time_exceeded or user_stop
            failed |= error_in_routine
            self.finished = run_finished
            self.failed   = failed
        #end if
    #end def check_sim_status


    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files


    def app_command(self):
        return self.app_name+' -input '+self.infile
    #end def app_command
#end class Pwscf



def generate_pwscf(**kwargs):
    sim_args,inp_args = Pwscf.separate_inputs(kwargs)

    if not 'input' in sim_args:
        input_type = inp_args.delete_optional('input_type','generic')
        sim_args.input = generate_pwscf_input(input_type,**inp_args)
    #end if
    pwscf = Pwscf(**sim_args)

    return pwscf
#end def generate_pwscf
