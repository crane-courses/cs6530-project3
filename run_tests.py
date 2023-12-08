import subprocess

start_size = 3
end_size = 13
for i in range(start_size, end_size + 1):
	exp_size = 2**i
	print('exp size: ' + str(exp_size))
	subprocess.call(['./test_bskip', str(exp_size), str('100000000')])
