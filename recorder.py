import sounddevice as sd
from scipy.io.wavfile import write


fs = 44100
sec = 5

i = 0
while True:
    answer = input("Record? (Y/N): ")

    if answer.upper() == 'Y':
        print("Recording...")
        record = sd.rec(int(fs * sec), samplerate=fs, channels=2, device=2)
        sd.wait()

        write("recording" + str(i) + ".wav", fs, record)
        print("Recording saved!")

        i += 1
    elif answer.upper() == 'N':
        break
    else:
        print("Your answer is not valid. Please try again.")