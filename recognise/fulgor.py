import os
import re
import struct
import sys

from collections import Counter


class Fulgor:

	def __init__(self, dbfile):

		self.cogs = {}

		with open(dbfile, "rb") as _in:
			n_cogs, *_ = struct.unpack('i', _in.read(4))
			for i in range(n_cogs):
				cog, n_colors = struct.unpack('ii', _in.read(8))
				# colors = struct.unpack(f'{n_colors}i', _in.read(4 * n_colors))
				# self.cogs[cog] = [colors, []]
				self.cogs[cog] = []

				n_csets, *_ = struct.unpack('i', _in.read(4))
				for j in range(n_csets):
					n, *_ = struct.unpack('i', _in.read(4))
					cset = struct.unpack(f'{n}i', _in.read(4 * n))
					# self.cogs[cog][1].append(cset)
					self.cogs[cog].append(cset)


	def eval_kcon(self, genome_id, fn):

		def parse_csets(s):
			for match in re.finditer(r'([0-9]+) ([0-9]+) ([0-9]+)', s):
				yield tuple(int(g) for g in match.groups())

		with open(fn) as _in:
			for line in _in:
				# 456320.20131.ABHB01000001_1219.COG0012  13      (0 33 249405)   (33 1 431313)   (34 53 249405)  (87 5 425771)   (92 349 249405) (441 8 335660)  (449 241 249405)        (690 17 249425) (707 232 249405)        (939 2 249425)  (941 175 249405)        (1116 5 610595) (1121 34 249405)
				p = line.find("\t")
				seqid, line = line[:p], line[p + 1:]
				cog = int(seqid[seqid.rfind(".") + 4:])
				line = line[line.find("\t") + 1:]
				#print(line)
				#seqid, _, *csets = line.rstrip().split("\t")
				c = Counter()

				n_kmers = 0
				for start, n, cset in parse_csets(line):
					# members = self.cogs[cog][1][cset]
					members = self.cogs[cog][cset]
					# print(start, n, cset, members)
					c.update(members * n)
					n_kmers += n

				# print(n_kmers, c)
				for k, v in c.items():
					fr = v/n_kmers
					if fr > 0.5:
						print(genome_id, seqid, k, v, fr)


def main():
	fulgor = Fulgor(sys.argv[1])

	# dirpath, dirs, files = os.walk(sys.argv[2])
	# for d in dirs:


	for d in sys.stdin:
		d = d.rstrip()
		f = f"{sys.argv[2]}/{d}/{d}.fulgor.out"
		if os.path.isfile(f):
			fulgor.eval_kcon(d, f)


	# fulgor.eval_kcon(sys.argv[2]) #, int(os.path.basename(f).replace("COG", "").replace(".out", "")))






if __name__ == "__main__":
	main()