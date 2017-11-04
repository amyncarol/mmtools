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
		if self.layer == 1:
			for i in os.listdir(basedir):
				folder = os.path.join(wd, i)
				dst_file = os.path.join(folder, 'submit.sh')
				copyfile(source_file, dst_file)
				self.submit_command()

		if self.layer == 2:
			for i in os.listdir(basedir):
				folder = os.path.join(wd, i)
				for j in os.listdir(folder):
					subfolder = os.path.join(folder, j)
					dst_file = os.path.join(subfolder, 'submit.sh')
					copyfile(source_file, dst_file)
					self.submit_command()









