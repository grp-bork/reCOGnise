import glob
import os
import struct
import sys

def read_colors(f):
	with open(f, "rt") as _in:
		cog = None
		for line in _in.readlines():
			_, cog_str, speci_file = line[line.find(" ") + 1:-1].split("/")
			if cog is None:
				cog = int(cog_str[3:])
				yield cog
			color = os.path.splitext(speci_file)[0]
			yield int(color[color.rfind("_") + 1:])

def read_colorsets(f, colors):
	with open(f, "rt") as _in:
		for line in _in:
			setid, _, *members = line.rstrip().split(" ")
			yield tuple(colors[int(m)] for m in members)



# def read_colorsets(fn, colors):
#     with open(fn, "rt") as _in:
#         for line in _in:
#             setid, _, *members = line.rstrip().split(" ")
#             # color_set_id=0 size=2 1005 1678
#             # yield int(setid.split("=")[1]), tuple(colors[int(m)] for m in members)
#             yield tuple(colors[int(m)] for m in members)

# color_set_id=0 size=2 31767 39655
# color_set_id=1 size=5 1246 7026 8800 8923 30708
# color_set_id=2 size=2 821 879
# color_set_id=3 size=2 582 22509
# color_set_id=4 size=2 13335 22370


def main():

	raw_db = sys.argv[1]

	color_files = glob.glob(f"{raw_db}/colors/COG*.colors.txt")

	with open(sys.argv[2], "wb") as _out:

		_out.write(struct.pack('i', len(color_files)))		

		for f in color_files:
			cog, *colors = read_colors(f)
			_out.write(struct.pack('ii', cog, len(colors)))
			# _out.write(struct.pack(f'{len(colors)}i', *colors))

			ff = f"{raw_db}/colorsets/COG{cog:04d}.color_sets.txt"
			colorsets = list(read_colorsets(ff, colors))
			_out.write(struct.pack('i', len(colorsets)))
			for cset in colorsets:
				_out.write(struct.pack(f'{len(cset) + 1}i', len(cset), *cset))


	
		



	# for f in COG0*.colors.txt; do color=$(basename $f .colors.txt); cut -f 2 $f | cut -f 3 -d / | cut -f 3 -d _ | sed "s/0*\([0-9]\+\)\.fa/\1/" | tr "\n" "," | sed "s/,$/\n/;s/^/"$color"\t/"; done > cog_speci_colors.txt

	# 0	cog_speci/COG0016/specI_v4_00000.fa
	# 1	cog_speci/COG0016/specI_v4_00001.fa
	# 2	cog_speci/COG0016/specI_v4_00002.fa
	# 3	cog_speci/COG0016/specI_v4_00003.fa
	# 4	cog_speci/COG0016/specI_v4_00004.fa

if __name__ == "__main__":
	main()