README

Algorytm mojego rozwiązania wygląda następująco:

-Dla każdej sekwencji tworzymy sekwencję komplementarną.

-Dla każdej sekwencji (zarówno teh wczytanej z pliku, jak i komplementarnej) 
tworzę jej klucz poprzez odpowiednie, podane w treści zadania operacje na 3-gramach. 
Tworzenie kluczy zrównoleglam dla każdej sekwencji.

-Sekwencje sortuję po ich kluczach.

-Korzystając z programowania równoległego tworzę graf powiązań sekwencji. Równocześnie uruchamiam  
obliczenia dla 40000 bloków. W każdym z nich równolegle odpalam algorytm needlemana-wunscha 
dla 128 par sekwencji. W obrębie bloku liczę wagę podobieństwa pomiędzy 
jedną z sekwencji a 128-moma kolejnymi na liście. Jako że pierwsza, analizowana w bloku, sekwencja 
wczytywana jest 128 razy, skopiowałam ją do pamięci współdzielonej (shared memory). 