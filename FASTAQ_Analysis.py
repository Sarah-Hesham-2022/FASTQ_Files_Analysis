import math
import matplotlib.pyplot as plot

file = open("File.fastq")
ids = []
seqs = []
scores = []
line = file.readline().rstrip()
ids.append(line)

while(line):
    line = file.readline().rstrip()
    seqs.append(line)
    line = file.readline().rstrip()
    line = file.readline().rstrip()
    scores.append(line)
    line = file.readline().rstrip()
    if(line):
        ids.append(line)

#print(len(seqs))
#print(len(ids))
#print(len(scores))
#print(seqs[2])
#print(ids[3])
#print(scores[4])

def Phred33LetterToErrorProbability(letterScore):

    ProbabilityOfError = int(letterScore) - 33
    ProbabilityOfError /= -10
    ProbabilityOfError = 10**ProbabilityOfError
    return ProbabilityOfError

def ErrorProbabilityToPhred33Letter(ProbabilityOfError):

    letterScore = (-10 * math.log10(ProbabilityOfError)) + 33
    return int(letterScore)

#print(Phred33LetterToErrorProbability(78))
#print(ErrorProbabilityToPhred33Letter(Phred33LetterToErrorProbability(78)))

#Exercise 1:
#Write a python script that reads and parses a FASTQ file by:
#Splitting the FASTQ file data into reads such that the reads are
#stored in a dictionary where the key is the sequence ID and the
#value is a list containing the sequence and the quality scores

myDict = {}

for i in range(len(ids)):
    myDict[ids[i]] = [seqs[i],scores[i]]

#print(myDict)
#print(len(myDict))

#Exercise 2:
#Write a function QtoPhred33 that takes Q as an input
#and returns the ASCII-encoded quality.
#Hint: chr() built-in function converts an integer to a character according to the ASCII-Table.
def QtoPhred33(Q):
    return chr(int(Q)+33)

#print(QtoPhred33(45.2))

#Write a function Phred33toQ that takes ASCII-encoded
#quality as an input and returns Q.
#Hint: ord() built-in function converts a character to an integer according to the ASCII-Table.
def Phred33toQ(letter):
    return (ord(letter)-33)

#print(Phred33toQ("N"))

#Exercise 3:
#In the previously parsed FASTQ File:
#For each read, calculate and print the quality score of each base
#(using Phred+33 scheme) and the accuracy of each base call.
#For each read, calculate the average quality and average accuracy.
readsBasesQScores = []
readsBasesAccuracy = []

for i in scores:
    temp1 = []
    temp2 = []

    for j in range(len(i)):
        temp1.append(Phred33toQ(i[j]))
        temp2.append(1-Phred33LetterToErrorProbability(ord(i[j])))

    readsBasesQScores.append(temp1)
    readsBasesAccuracy.append(temp2)

#print(len(readsBasesQScores))
#print(len(readsBasesAccuracy))
#print(len(readsBasesQScores[2]))
#print(len(readsBasesAccuracy[2]))
#print(readsBasesQScores[4])
#print(readsBasesAccuracy[4])

avgReadQScores = []
avgReadAccuracies = []

for i,m in zip(readsBasesAccuracy, readsBasesQScores):
    sum1 = 0
    sum2  = 0

    for j,n in zip(range(len(i)), range(len(m))):
        sum1 = sum1 + i[j]
        sum2 = sum2 + m[n]

    avgReadAccuracies.append(sum1/len(i))
    avgReadQScores.append(int(sum2/len(m)))

#print(avgReadAccuracies)
#print(avgReadQScores)
#print(len(avgReadAccuracies))
#print(len(avgReadQScores))

#Let’s build and plot a histogram of base qualities read from the given FASTQ file.
hist = []
uniqueScores = set(avgReadQScores)
uniqueScores = list(uniqueScores)

for i in uniqueScores:
    hist.append(avgReadQScores.count(i))

#print(uniqueScores)
#print(hist)
#plot.xlabel('Quality Score')
#plot.ylabel('Frequency')
#plot.title('Base Quality Histogram')
#plot.bar(uniqueScores, hist)
#plot.show()
#plot.plot(uniqueScores, hist)
#plot.show()

#Exercise 1
#Write a python script that reads the given FASTQ file and:
#Calculates the “Per base sequence quality” (mean and median only)
#Hint: Use sorted() built-in function to sort a list and easily get the median.
#Draws a plot (any kind of plot you prefer) for this metric
#Open the same FASTQ file in FastQC and compare the results of “Per base sequence quality” to your output

avgPerBaseScore = []
medianPerBaseScores = []
seqsLengths = []
maxLen = 0

for i in seqs:
    seqsLengths.append(len(i))

maxLen = sorted(seqsLengths)[-1]

for i in range(maxLen):
    sum = 0
    count = 0
    for j in range(len(readsBasesQScores)):
        try:
            sum = sum + readsBasesQScores[j][i]
            count = count + 1
        except IndexError:
            continue
    avgPerBaseScore.append(int(sum/count))

#print(avgPerBaseScore)
#print(len(avgPerBaseScore))

for i in range(maxLen):
    temp = []
    for j in range(len(readsBasesQScores)):
        try:
            temp.append(readsBasesQScores[j][i])
        except IndexError:
            continue
    temp = sorted(temp)
    median = temp[int(len(temp)/2)]
    if(len(temp)%2 == 0):
         median = (temp[int(len(temp)/2 - 1)] + temp[int(len(temp)/2)])/2
    medianPerBaseScores.append(int(median))

#print(medianPerBaseScores)
#print(len(medianPerBaseScores))

#➔ In the previously used FASTQ file, write a python script that:
#Calculates the “Per base sequence content”
#Draws a plot (any kind of plot you prefer) for this metric
#Open the same FASTQ file in FastQC and compare the results of “Per base sequence content” to your output.

content_score = []
a= []
c= []
g= []
t= []
positions = []

for i in range(maxLen):
    temp=[]
    for j in range(len(seqs)):
        try:
            temp.append(seqs[j][i])
        except IndexError:
            continue
    positions.append(i+1)
    a.append(int(temp.count("A")/len(temp)*100))
    c.append(int(temp.count("C")/len(temp)*100))
    g.append(int(temp.count("G")/len(temp)*100))
    t.append(int(temp.count("T")/len(temp)*100))
    content_score.append([int(temp.count("A")/len(temp)*100) ,int(temp.count("C")/len(temp)*100) ,int(temp.count("G")/len(temp)*100) ,int(temp.count("T")/len(temp)*100)])

#print(len(content_score))
#print(len(content_score[2]))
#print(content_score[5])
#plot.xlabel('Position in Read')
#plot.ylabel('Percentage of Nucleotides')
#plot.title('Per Base Sequence Content')
#plot.plot(positions,a,color='green')
#plot.plot(positions,c,color='blue')
#plot.plot(positions,g,color='black')
#plot.plot(positions,t,color='red')
#plot.ylim(0,100)
#plot.xlim(1,151)
#plot.show()

#Per Sequence GC Content

gc = []
for i in seqs:
    gc.append(int((i.count("G")+i.count("C"))/len(i)*100))
 
#print(len(gc))
#print(gc)

counts = []
uniqueCounts = set(gc)
uniqueCounts = list(uniqueCounts)
for i in uniqueCounts:
    counts.append(gc.count(i))

#print(uniqueCounts)
#print(counts)

#plot.xlabel('Percentages of Per Sequence GC Content')
#plot.ylabel('Count')
#plot.title('Per Sequence GC Content')
#plot.plot(uniqueCounts,counts,color='red')
#plot.show()

#Per Base N Content

nbases = []
for i in range(maxLen):
    count = 0
    sum = 0
    for j in range(len(seqs)):
        try:
            count = count +1
            if(seqs[j][i] == "N" or seqs[j][i]=="n"):
                sum = sum +1
        except IndexError:
            continue
    nbases.append(int(sum/count*100))

#print(len(nbases))
#print(nbases)

#plot.xlabel('Percentages of Per Sequence N Content')
#plot.ylabel('Position in Read')
#plot.title('Per Base N Content')
#plot.plot(positions,nbases,color='red')
#plot.ylim(0,100)
#plot.show()

#Sequence Length Distribution

uniqueLen = set(seqsLengths)
uniqueLen = list(uniqueLen)
countLen = []

for i in uniqueLen:
    countLen.append(seqsLengths.count(i))

#print(countLen)
#print(uniqueLen)
#print(len(seqsLengths))

#plot.ylabel('Count')
#plot.xlabel('Position in Read')
#plot.title('Sequence Length Distribution')
#plot.plot(uniqueLen,countLen,color='red')
#plot.bar(uniqueLen,countLen,color='red')
#plot.xlim(1,2*maxLen+1)
#plot.ylim(0,2*len(seqs)+1)
#plot.show()

#Sequence Duplication Levels

repeatDict = {}
uniqueSeqs = set(seqs)
uniqueSeqs = list(uniqueSeqs)

duplications = []
for i in uniqueSeqs:
    duplications.append(seqs.count(i))
    repeatDict[i] = duplications[-1]

#print(duplications)
#print(len(uniqueSeqs))
#print(len(duplications))
#print(len(seqs))

uniqueDuplications = set(duplications)
uniqueDuplications = list(uniqueDuplications)
duplicationsCount = []
total = 0

for i in uniqueDuplications:
    duplicationsCount.append(duplications.count(i))
    total = total + duplicationsCount[-1]

#print(uniqueDuplications)
#print(duplicationsCount)

percentagesDuplicationsCount = []

for i in duplicationsCount:
    percentagesDuplicationsCount.append(int(i/total*100))

#print(percentagesDuplicationsCount)

#plot.xlabel('Duplications Levels')
#plot.ylabel('Percentage of Duplications')
#plot.title('Sequence Duplication Levels')
#plot.plot(uniqueDuplications,percentagesDuplicationsCount,color='red')
#plot.ylim(0,100)
#plot.show()

#Overrepresented Sequences, This module shows the list of sequences which appear
#more than expected in the file.
#A sequence is considered overrepresented if it accounts for ≥ 0.1% of the total reads.

#print(repeatDict)
#print(len(repeatDict))

overRepresentedSeqs = []

for i in repeatDict:
    rate = round((repeatDict[i]/len(seqs)),3)*100
    if( rate > 0.1):
        overRepresentedSeqs.append(i)

#print(overRepresentedSeqs)


