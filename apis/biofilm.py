import sys, getopt

def main(argv):
	feature = 0
	predict = 0
	input_file_path = ""
	output_file_path = ""
	test_file_path = ""
	SVM_joblib_file_path = ""
	str_help = "biofilm USAGE:\n  biofilm.py -f <feature number> -p <perform prediction> -t <test file path for prediction> -i <input file path> -o <output file path>\n" +\
	"\n Please select features from the list below: \n  1- AAC\n  2- DPC\n  3- CTD\n"+\
	"\n If you want to perform prediction set the value 1 for -p: \n  -p 1"
	try:
		opts, args = getopt.getopt(argv, "hf:i:o:p:t:j:", ["feature=", "predict=", "input=", "output=", "test=", "joblib="])
	except getopt.GetoptError:
		print(str_help)
		sys.exit()
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print(str_help)
			sys.exit()
		if opt in ("-f", "--feature"):
			try:
				feature = int(arg)
			except Exception as e:
				print(str_help + "\n   Error: -f should be an Integer")
				sys.exit()
		if opt in ("-p", "--predict"):
			try:
				predict = int(arg)
			except Exception as e:
				print(str_help + "\n   Error: -p should be an Integer")
				sys.exit()
		if opt in ("-i", "--input"):
			try:
				input_file_path = arg
			except Exception as e:
				print(str_help + "\n   Error: -i should be a String")
				sys.exit()
		if opt in ("-o", "--output"):
			try:
				output_file_path = arg
			except Exception as e:
				print(str_help + "\n   Error: -o should be a String")
				sys.exit()
		if opt in ("-t", "--test"):
			try:
				test_file_path = arg
			except Exception as e:
				print(str_help + "\n   Error: -t should be a String")
				sys.exit()
		if opt in ("-j", "--joblib"):
			try:
				SVM_joblib_file_path = arg
			except Exception as e:
				print(str_help + "\n   Error: -j should be a String")
				sys.exit()

	#for Feature extraction
	if feature == 1:
		import AAC1
		AAC1.CalculateAAC4All(input_file_path, output_file_path)

	if feature == 2:
		import DPC
		DPC.CalculateDPC4All(input_file_path, output_file_path)

	if feature == 3:
		import CTD1
		CTD1.CalculateCTD4All(input_file_path, output_file_path)

	if predict == 1:
		import prediction
		prediction.perform_prediction(SVM_joblib_file_path, test_file_path, output_file_path)





if __name__ == "__main__":
	main(sys.argv[1:])













#
