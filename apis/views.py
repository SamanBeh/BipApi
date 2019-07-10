from django.shortcuts import render
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
from apis.models import Car
import json
from prediction import perform_prediction

def index(request):
    response = json.dumps([{}])
    return HttpResponse(response, content_type='text/json')

def get_car(request, car_name):
    if request.method == 'GET':
        try:
            car = Car.objects.get(name=car_name)
            response = json.dumps([{'Car': car.name, 'Top Speed': car.top_speed}])
        except:
            response = json.dumps([{"Error": "No car with that name"}])
    return HttpResponse(response, content_type='text/json')

@csrf_exempt
def add_car(request):
    if request.method == 'POST':
        payload = json.loads(request.body)
        car_name = payload['car_name']
        top_speed = payload['top_speed']
        car = Car(name=car_name, top_speed=top_speed)
        try:
            car.save()
            response = json.dumps([{'Success': 'Car added successfully!'}])
        except:
            response = json.dumps([{'Error': 'Car could not be added!'}])
    return HttpResponse(response, content_type='text/json')

@csrf_exempt
def add_car(request):
    if request.method == 'POST':
        payload = json.loads(request.body)
        lines = payload['lines']

        from joblib import load
        import pandas as pd
        import json
        # testdf = pd.DataFrame([[1,2],[3,4]], columns=["a", "b"])
        # testdf.to_csv("D:\\vhosts\\cbb1.ut.ac.ir\\httpdocs\\UploadedFiles\\testdf.csv", sep=",", index=False)

        # lines = json.loads(test_file_lines)
        df = pd.DataFrame(columns = list(lines[0].split(",")))
        for line, ind in zip(lines[1:], range(len(lines)-1)):
            df.loc[ind] = line

        model = load('SVMModel.joblib') # 'SVMModel.joblib'
        # X_test = pd.read_csv(test_file_path)
        X_test = df
        predictions = model.predict(X_test.iloc[:,1:])
        df = pd.DataFrame(columns=["Peptide sequence", "Biofilm inhobitor"])
        seqs = list(X_test.seq)
        for indx in range(len(predictions)):
            df = df.append({"Peptide sequence": seqs[indx], "Biofilm inhobitor": str(predictions[indx]).replace("0", "non BIP").replace("1", "BIP")}, ignore_index=True)
        df.to_csv(output_file_path, sep=',', index=False)
        # print(X_test)
        # print(predictions)
        # print(len(predictions))















#
