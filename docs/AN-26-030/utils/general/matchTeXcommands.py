import re
import os
import string
import glob


def get_cmd_matchstr(cmd):
	return '('+cmd.replace('\\','\\\\')+')'+'[^a-zA-Z]' # this would miss a macro with a non-alpha character (e.g. @)


def find_all_deps(cmddefs,cmd,used):
	cmd_match = '(\\\\[a-zA-Z]*)'
	m = re.findall(cmd_match,cmddefs[cmd])
	if m:
		for ctmp in m:
			if not ctmp or ctmp==cmd:
				continue
			if ctmp in cmddefs.keys():
				if ctmp not in used:
					used.append(ctmp)
				find_all_deps(cmddefs,ctmp,used)


def main():
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument("texfiles", help="names of TeX files to check", type=str, nargs="+")
	parser.add_argument("--verbose", help="verbose printout", action='store_true', required=False, default=False)
	parser.add_argument("--stop", type=str, default="cmsNoteHeader", help='stop recognizing commands after first line matching this string (to disable, use --stop "")')
	args = parser.parse_args()
	texfiles = args.texfiles
	verbose = args.verbose

	match_prefix = '\\\\(?:(?:re)*new|provide)command\s*\{?\s*' # might consider looking for \def as well
	match_suffix = '\s*\}?\s*(?:\[.*\])?\s*\{(.*)\}'
	# That means you should NEVER iterate the tex files line-by-line!
	p = re.compile(('(' + match_prefix + '(.*?)' + match_suffix + ')'))
	cmdorder = []
	cmds = dict()
	cmddefs = dict()
	cmdfulldefs = dict()
	# first pass: find command definitions
	for texfile in texfiles:
		with open(texfile) as tex:
			body = tex.read()
			# restrict first pass if requested
			body1 = body
			if len(args.stop)>0:
				stop_index = body.find(args.stop)
				if stop_index<0:
					print("Warning: requested stop string {} not found, will search entire file".format(args.stop))
				else:
					body1 = body[:stop_index]
			m = p.findall(body1)
			for c in m:
				if c[0] is not None and c[0]!='' and c[1] is not None and c[1]!='' and c[1] not in cmdorder:
					cmdorder.append(c[1])
					cmddefs[c[1]] = c[2]
					cmdfulldefs[c[1]] = c[0]

	cmdorder.reverse()

	# The total count of a command equals
	# D=1 for its declaration
	# + N number of times it is used in nested commands
	# + M number of times it is used in the main text.
	# When we make the second pass later, the count will include all of the above.
	# What we are actually interested is only M though,
	# so we subtract (1+N) here.
	for cmd in cmdorder:
		neg_count = 1 # D
		cmd_match = get_cmd_matchstr(cmd)
		for cc in cmddefs.keys():
			m = re.findall(cmd_match,cmddefs[cc])
			if m:
				neg_count += len(m) # N = sum_file (N_file)
		cmds[cmd] = -neg_count # Subtract!

	# second pass: find usage
	used = []
	for cmd in cmdorder:
		for texfile in texfiles:
			with open(texfile) as tex:
				body = tex.read()
				cmd_match = get_cmd_matchstr(cmd)
				m = re.findall(cmd_match,body)
				if m:
					cmds[cmd] += len(m) # Simply add!
		if cmds[cmd]<0:
			raise RuntimeError('{} count is {}<0.'.format(cmd,cmds[cmd]))
		elif cmds[cmd]>0:
			used.append(cmd)

   # This is where it becomes useful to have reversed cmdorder.
	for i in range(len(cmdorder)):
		cmd = cmdorder[i]
		if cmd not in used:
			continue
		find_all_deps(cmddefs,cmd,used)

	cmdorder.reverse()

	if verbose:
		used_in_text = [cmd for cmd in cmdorder if cmds[cmd]>0]
		used_onlynested = [cmd for cmd in cmdorder if cmds[cmd]==0 and cmd in used]
		print("----------------\nUser-defined commands found in the main TeX files:\n----------------")
		for u in [cmd for cmd in cmdorder if cmd in used_in_text]:
			print('{:>2} occurrences of {}'.format(cmds[u], u))
		print("----------------\nCommands only nested in other user-defined commands:\n----------------")
		for  u in [cmd for cmd in cmdorder if cmd in used_onlynested]:
			print(u)
	print("----------------\nCommands unused anywhere:\n----------------")
	for u in [cmd for cmd in cmdorder if cmd not in used]:
		print(u)
	if verbose:
		print("----------------\nFinal list of commands to keep:\n----------------")
		for  u in [cmd for cmd in cmdorder if cmd in used]:
			print(u)
	if verbose:
		print("----------------\nFinal list of commands to include in LaTeX form:\n----------------")
		for  u in [cmd for cmd in cmdorder if cmd in used]:
			print(cmdfulldefs[u])


if __name__ == "__main__":
	main()
