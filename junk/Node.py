def addNodes(fakegene):
	N7 = "N7"
	N6 = "N6"
	N5 = "N5"
	N4 = "N4"
	N3 = "N3"
	N2 = "N2"
	N1 = "N1"
	N0 = "N0"

	List = []
	for i in fakegene:
		if i == "gene":
			List.append(1)
		else:
			List.append(0)

	# Look at N7, N6
	if List[0] + List[1] > 1:
		N7 = N7
		if List[2] == 1:
			N6 = N6
		else:
			N3 = N3 + N6
			N6 = ""
	else:
		N6 = N6 + N7
		N7 = ""
		if List[0] + List[1] + List[2] > 1:
			N6 = N6
		else:
			N3 = N3 + N6
			N6 = ""

	# Look at N5, N4
	if List[3] + List[4] > 1:
		N5 = N5
		if List[5] == 1:
			N4 = N4
		else:
			N3 = N3 + N4
			N4 = ""
	else:
		N4 = N4 + N5
		N5 = ""
		if List[3] + List[4] + List[5] > 1:
			N4 = N4
		else:
			N3 = N3 + N4
			N4 = ""

	# Look at N3
	if List[0] + List[1] + List[2] > 0 and List[3] + List[4] + List[5] > 0:
		N3 = N3
	else:
		N1 = N1 + N3
		N3 = ""

	# Look at N2
	if List[6] + List[7] > 1:
		N2 = N2
	else:
		N1 = N1 + N2
		N2 = ""

	# Look at N1
	if List[0] + List[1] + List[2] > 0:
		N3_1 = 1
	else:
		N3_1 = 0
	if List[3] + List[4] + List[5] > 0:
		N3_2 = 1
	else:
		N3_2 = 0
	if List[6] + List[7] > 0:
		N2_1 = 1
	else:
		N2_1 = 0
	if N3_1 + N3_2 + N2_1 > 1:
		N1 = N1
	else:
		N0 = N0 + N1
		N1 = ""



	nodeList = [N7, N6, N5, N4, N3, N2, N1, N0]
	output = []
	for i in nodeList:
		if i != "":
			output.append(i)
	return output


if __name__ == '__main__':
	fakegene = ['', '', 'gene', '', 'gene', '', '', 'gene']

	print(addNodes(fakegene))

