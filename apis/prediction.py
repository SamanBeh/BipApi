def dataframe_to_json(df):
    import json
    dct = {}
    for cn in list(df):
        dct[cn] = list(df[cn])

    return json.dumps(dct)


def perform_prediction(SVM_joblib_file_path, test_file_path, output_file_path):
    #for predicting
    #TODO: extract features | merge them | then perform prediction
    from joblib import load
    with open(output_file_path, "a") as f:
        f.write("YAAAYYY joblib")
    import pandas as pd
    # testdf = pd.DataFrame([[1,2],[3,4]], columns=["a", "b"])
    # testdf.to_csv("D:\\vhosts\\cbb1.ut.ac.ir\\httpdocs\\UploadedFiles\\testdf.csv", sep=",", index=False)
    with open(output_file_path, "w") as f:
        f.write("YAAAYYY pandas imported ")
    model = load(SVM_joblib_file_path) # 'SVMModel.joblib'
    with open(output_file_path, "a") as f:
        f.write("YAAAYYY SVMModel.joblib loaded successfully")
    X_test = pd.read_csv(test_file_path)
    predictions = model.predict(X_test.iloc[:,1:])
    df = pd.DataFrame(columns=["Peptide sequence", "Biofilm inhobitor"])
    seqs = list(X_test.seq)
    for indx in range(len(predictions)):
        df = df.append({"Peptide sequence": seqs[indx], "Biofilm inhobitor": str(predictions[indx]).replace("0", "non BIP").replace("1", "BIP")}, ignore_index=True)
    df.to_csv(output_file_path, sep=',', index=False)
    # print(X_test)
    # print(predictions)
    # print(len(predictions))
