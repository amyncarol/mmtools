import mmtools.data as data
import os
from shutil import copyfile
import subprocess
class JobRunner(object):
	"""
	generates the submit.sh files for thor, fenrir and stampede
	"""
	def __init__(self, nodes, cores, vasp_exe_type, cluster):
		"""
		vasp_exe_type: can be std, gam, ncl
		cluster: can be thor, fenrir, or stampede
		"""
		self.nodes = nodes
		self.cores = cores
		self.vasp_exe_type = vasp_exe_type
		self.cluster = cluster

	def submit_file(self):
		"""
		generate a single submit.sh file
		"""
		path = os.path.dirname(data.__file__)
		wd = os.getcwd()
		lines = []
		if self.cluster == 'thor':
			with open(os.path.join(path, 'submit_thor.sh'),'r') as f:
				for line in f:
					if line == '#PBS -l nodes=1:ppn=16\n':
						line = '#PBS -l nodes=1:ppn='+str(self.cores)+'\n'
					if line == "MYMPIPROG=\"${HOME}/vtst/bin/vasp_std\"\n":
						line = "MYMPIPROG=\"${HOME}/vtst/bin/" + "vasp_" + self.vasp_exe_type +"\"\n"
					lines.append(line)
		with open(os.path.join(wd, 'submit.sh'),'w') as f:
			for line in lines:
				f.write(line)
	##TODO: for fenir and stampede

	def submit_command(self):
		if self.cluster == 'thor' or 'fenrir':
			command = "qsub submit.sh"
		elif self.cluster == 'stampede':
			command = "sbatch submit.sh"
		p = subprocess.Popen(command.split())

	def batch_runner(self, layer):
		"""
		layer: can be 1 or 2, 2 means the running directory is 2 layer down current directory
		"""
		self.submit_file()
		wd = os.getwd()
		source_file = os.path.join(wd, 'submit.sh')
		
		for i in os.listdir(basedir):
			folder = os.path.join(wd, i)
			
			if self.layer == 2:
				for j in os.listdir(folder):
					subfolder = os.path.join(folder, j)
					dst_file = os.path.join(subfolder, 'submit.sh')
					copyfile(source_file, dst_file)
					self.submit_command()

			elif self.layer == 1:
				dst_file = os.path.join(folder, 'submit.sh')
				copyfile(source_file, dst_file)
				self.submit_command()


def bundle_submit_file(folder, path_of_custodian_file, n_jobs, queue='skx-normal', \
		n_node=4, n_tasks_per_node=24, walltime='48:00:00'):
	"""
	generate several submit scripts, each submit script can run multiple vasp jobs, this is a way to deal with
	the limit of number of jobs that can be submitted on some clusters.

	Args:
		folder: the parent folder of the vasp folders 
		n_jobs: the number of submit files to generate, should be less than the limit of jobs that can be submitted at a time.
		queue: after "#SBATCH -p"
		n_node: after "#SBATCH -N"
		n_tasks_per_node: after "#SBATCH --ntasks-per-node"
		walltime: after "#SBATCH -t"

	Writes:
		submit_0.sh
		submit_1.sh
		....

	Example submit file:
		#!/bin/bash

		#SBATCH -J myjob           
		#SBATCH -o myjob.o%j      
		#SBATCH -e myjob.e%j     
		#SBATCH -p skx-dev       
		#SBATCH -N 4               
		#SBATCH --ntasks-per-node 24             
		#SBATCH -t 00:30:00        
		###SBATCH -A myproject       

		ibrun /home1/05018/tg843171/vasp.5.4.4_vtst/bin/vasp_std
	"""
	subfolders = [subfolder for subfolder in os.listdir(folder) if os.path.isdir(os.path.join(folder, subfolder))]
	n_batch = len(subfolders)//n_jobs + 1 * (len(subfolders)%n_jobs != 0)
	for i in range(n_jobs):
		with open(os.path.join(folder, 'submit_' + str(i) + '.sh'), 'w') as f:
			f.write('#!/bin/bash\n')
			f.write('#SBATCH -J myjob\n')
			f.write('#SBATCH -o myjob.o%j\n')      
			f.write('#SBATCH -e myjob.e%j\n')    
			f.write('#SBATCH -p ' + queue+'\n')      
			f.write('#SBATCH -N ' + str(n_node)+'\n')             
			f.write('#SBATCH --ntasks-per-node ' + str(n_tasks_per_node)+'\n')  
			f.write('#SBATCH -t '+walltime+'\n')        
			f.write('###SBATCH -A myproject\n') 
			for j in range(i*n_batch, min((i+1)*n_batch, len(subfolders))):
				f.write('cd ' + subfolders[j]+'\n')
				f.write('cp ' + path_of_custodian_file + ' ./custodian_job.py\n')
				f.write('python custodian_job.py\n')
				f.write('cd ..\n')

if __name__=='__main__':
	bundle_submit_file('/Users/yao/Google Drive/mmtools/workflow/test_data', \
		'/Users/yao/Google Drive/mmtools/workflow/test_data/custodian_job.py', 2)






