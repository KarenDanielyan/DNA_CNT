#Imports
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import HashingVectorizer
import re
import pandas as pd

# Reading the file
set = pd.read_csv(r'C:\Users\User\Desktop\DNA_CNT\Data\seq_TC_12mer_2strun.csv')
set.columns = ['Sequence', 'Class']
seq =  set['Sequence']

# Data preprocessing
# Removing the non-nucleotide characters
seq = seq.apply(lambda x: re.sub('[^ACGT]', '', x))
# Removing the empty sequences
seq = seq.apply(lambda x: x.strip())

# Count vectorization 
# Creating the vectorizer
vectorizer = CountVectorizer(analyzer='char', ngram_range=(1,3))
# Fitting the vectorizer
vectorizer.fit(seq)

# Hash Vectorization 
# Creating the vectorizer
hash_vectorizer = HashingVectorizer(analyzer='char', ngram_range=(1,3))
# Fitting the vectorizer
hash_vectorizer.fit(seq)

# TF-IDF Vectorization
# Creating the vectorizer
term_freq_vectorizer = TfidfVectorizer(analyzer='char', ngram_range=(1,3))
# Fitting the vectorizer
term_freq_vectorizer.fit(seq)
