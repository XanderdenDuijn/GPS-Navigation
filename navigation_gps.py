with open(fn_monuments,'r') as json_monuments, open(fn_tweets,'r') as json_tweets:
        #receiving the dicts
        data_monuments = json.load(json_monuments)
        data_tweets = json.load(json_tweets)
