import pickle

with open("GPICD59/TEST", "rb") as fp:
    b = pickle.load(fp)

print(b)