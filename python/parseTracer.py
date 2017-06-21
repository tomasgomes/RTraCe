from sys import argv

infile1 = argv[1]
infile2 = argv[2]
outfile = argv[3]

# parse recombinants.txt
cells_dic = {}
with open(infile1) as f:
	for line in f:
		if "\t" in line and not "cell_name" in line:
			line = line.strip("\n").split("\t")
			cell = line[0]
			if not cell in cells_dic.keys():
				cells_dic[cell] = ["notAssigned", "No", "No", "No", "No",
				"No", "No", "No", "No",
				"No", "No", "No", "No",
				"No", "No", "No", "No"]
			if line[1].startswith("No seqs"):
				cells_dic[cell] = ["NoTCR", "NA", "NA", "NA", "NA",
				"NA", "NA", "NA", "NA",
				"No", "No", "No", "No",
				"No", "No", "No", "No"]

			elif line[1]=="A":
				seq = line[2]
				if line[3]=="True":
					prod = "Yes"
				else:
					prod = "No"
				if cells_dic[cell][1]=="No":
					cells_dic[cell][1] = seq
					cells_dic[cell][9] = prod
				else:
					cells_dic[cell][2] = seq
					cells_dic[cell][10] = prod

			elif line[1]=="B":
				seq = line[2]
				if line[3]=="True":
					prod = "Yes"
				else:
					prod = "No"
				if cells_dic[cell][3]=="No":
					cells_dic[cell][3] = seq
					cells_dic[cell][11] = prod
				else:
					cells_dic[cell][4] = seq
					cells_dic[cell][12] = prod

			elif line[1]=="G":
				seq = line[2]
				if line[3]=="True":
					prod = "Yes"
				else:
					prod = "No"
				if cells_dic[cell][5]=="No":
					cells_dic[cell][5] = seq
					cells_dic[cell][13] = prod
				else:
					cells_dic[cell][6] = seq
					cells_dic[cell][14] = prod

			elif line[1]=="D":
				seq = line[2]
				if line[3]=="True":
					prod = "Yes"
				else:
					prod = "No"
				if cells_dic[cell][7]=="No":
					cells_dic[cell][7] = seq
					cells_dic[cell][15] = prod
				else:
					cells_dic[cell][8] = seq
					cells_dic[cell][16] = prod

# add simplified clonotypes and other cell types
multi = False
iNKT = False
MAIT = False
clon = False
cl_n = 0
with open(infile2) as f:
	for line in f:
		if "Cells with" in line:
			multi = True
		elif "iNKT" in line:
			multi = False
			iNKT = True
		elif "MAIT" in line:
			iNKT = False
			MAIT = True
		elif "Clonotype" in line:
			iNKT = False
			MAIT = False
			clon = True

		if line.startswith("###"):

			cell = line.strip("### ").rstrip(" ###\n")
			if cell not in cells_dic.keys():
				cells_dic[cell] = ["notAssigned", "No", "No", "No", "No",
				"No", "No", "No", "No",
				"No", "No", "No", "No",
				"No", "No", "No", "No"]
			if multi:
				cells_dic[cell][0] = "Multi_recomb"
			elif iNKT:
				cells_dic[cell][0] = "iNKT"
			elif MAIT:
				cells_dic[cell][0] = "MAIT"

		if clon and "," in line:
			cl_n += 1
			line = line.strip("\n").split(", ")
			for cell in line:
				cells_dic[cell][0] = "cl" + str(cl_n)



# write output
header = "tracer_cell_name\tclonotype_simp\tA_1\tA_2\tB_1\tB_2\tG_1\tG_2\tD_1\tD_2" +\
"\tpA_1\tpA_2\tpB_1\tpB_2\tpG_1\tpG_2\tpD_1\tpD_2\n"
header = header.replace("\t", "\ttracer_")
with open(outfile, "w") as res:
	res.write(header)
	for key in cells_dic.keys():
		newline = key + "\t" + "\t".join(cells_dic[key])+"\n"
		res.write(newline)
