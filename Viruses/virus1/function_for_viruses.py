def seq_Nucleotide_composition(seq):
    N =  seq.count('N')
    n =  seq.count('n')
    a =  round(seq.count('a')/(len(seq)-N)*100,2)
    t =  round(seq.count('t')/(len(seq)-N)*100,2)
    c =  round(seq.count('c')/(len(seq)-N)*100,2)
    g =  round(seq.count('g')/(len(seq)-N)*100,2)
    A =  round(seq.count('A')/(len(seq)-N)*100,2)
    T =  round(seq.count('T')/(len(seq)-N)*100,2)
    C =  round(seq.count('C')/(len(seq)-N)*100,2)
    G =  round(seq.count('G')/(len(seq)-N)*100,2)


    #print('A',A,'T',T,'C',C,'G',G,'N',N)
    #print('a',a,'t',t,'c',c,'g',g,'n',n)
    #print('all nucleotide',a+t+c+g+A+T+C+G)
    return A,T,C,G



def calculat_portion_of_elements_overlap_sharing(sequence, gff_path,chr_id,element_names):
    # this function is for labling + and - and shered region in sequence
    # result labled sequence
    # 1234 for non maped
    # ATCG for +
    # atcg for -
    # ~!@# for shered

    len_seq = len(sequence)
    Element_number = 0
    upper_sequence = sequence.upper()#remove prore mask

    sequence_number_base = upper_sequence.replace('A','1')
    sequence_number_base = sequence_number_base.replace('T','2')
    sequence_number_base = sequence_number_base.replace('C','3')
    sequence_number_base = sequence_number_base.replace('G','4')
    sequence_number_base = list(sequence_number_base) #sequence with number 1234

    #gff = gzip.open(gff_path,mode='rt')
    gff = open(gff_path)

    #step 1 for plus
    element_counter = 0
    for annotation in gff:
        #print(annotation)
        if annotation[0] == '#':
            continue
        annotation = annotation.strip().split()
        start_point = 3
        end_point = 4
        orientation = 6
        element_point = 2
        if len(annotation)<5:
            continue
        if annotation[start_point] == '1':
            continue
        if annotation[3] == 'CDS':
            start_point += 1
            end_point += 1
            orientation += 1
            element_point +=1
        #print(annotation[2],annotation)
        for element in element_names:
            #print(element)
            repeat_control = []
            if element == annotation[element_point] and annotation[0] == chr_id:
                #print(element)
                element_counter += 1

                if int(annotation[end_point]) - int(annotation[start_point]) == 0:
                    '''print('len element was 0')
                    print(annotation)
                    print(len_seq)'''
                    continue
                Element_number += 1
                if annotation[orientation] == '+':
                    #print(1)
                    P_seq = sequence_number_base[int(annotation[start_point]) - 1:int(annotation[end_point])]
                    #print(start_point,end_point)
                    P_seq = ''.join(P_seq)
                    P_seq = P_seq.replace('1','A')
                    P_seq = P_seq.replace('2','T')
                    P_seq = P_seq.replace('3','C')
                    P_seq = P_seq.replace('4','G')
                    P_seq = P_seq.replace('a','~')
                    P_seq = P_seq.replace('t','!')
                    P_seq = P_seq.replace('c','@')
                    P_seq = P_seq.replace('g','#')
                    P_seq = list(P_seq)
                    #print('P',P_seq[:100])

                    sequence_number_base[int(annotation[start_point]) - 1:int(annotation[end_point])] = P_seq



                if annotation[orientation] == '-':
                    C_seq = sequence_number_base[int(annotation[start_point]) - 1:int(annotation[end_point])]
                    C_seq = ''.join(C_seq)
                    C_seq = C_seq.replace('1','a')
                    C_seq = C_seq.replace('2','t')
                    C_seq = C_seq.replace('3','c')
                    C_seq = C_seq.replace('4','g')
                    C_seq = C_seq.replace('A','~')
                    C_seq = C_seq.replace('T','!')
                    C_seq = C_seq.replace('C','@')
                    C_seq = C_seq.replace('G','#')
                    C_seq = list(C_seq)
                    #print('C',C_seq[0:100])
                    sequence_number_base[int(annotation[start_point]) - 1:int(annotation[end_point])] = C_seq


    if element_counter == 0 :
        #print('map error')
        return 'map error'
    sequence_number_base = ''.join(sequence_number_base)

    N = sequence_number_base.count('N')
    #print('N',N)
    nonmap_portion = (sequence_number_base.count('1') + sequence_number_base.count('2') + sequence_number_base.count('3') + sequence_number_base.count('4'))/len(sequence_number_base) * 100
    plus_portion = (sequence_number_base.count('A') + sequence_number_base.count('T') + sequence_number_base.count('C') + sequence_number_base.count('G'))/len(sequence_number_base) * 100
    minis_portion = (sequence_number_base.count('a') + sequence_number_base.count('t') + sequence_number_base.count('c') + sequence_number_base.count('g'))/len(sequence_number_base)*100
    share_portion = (sequence_number_base.count('~') + sequence_number_base.count('!') + sequence_number_base.count('@') + sequence_number_base.count('#'))/len(sequence_number_base) * 100
    '''print('Plus portion',plus_portion)
    print('Negative portion',minis_portion)
    print('share portion',share_portion)
    print('analysis eficincy (start len to finsh len)',len(sequence_number_base)/len(sequence)*100,'%')
    print('nonmap',nonmap_portion)
    print(element_names,Element_number)'''
    #print(sequence_number_base)
    #print(annotation)
    #print(chr_id)
    return sequence_number_base



def overlap_nucleotide_composition(sequence):
    #masked region nucleotide composition
    Len_seq = len(sequence)

    A = sequence.count('A')
    T = sequence.count('T')
    C = sequence.count('C')
    G = sequence.count('G')
    one_masked_count_P = A + T + C + G



    a = sequence.count('a')
    t = sequence.count('t')
    c = sequence.count('c')
    g = sequence.count('g')
    one_masked_count_C = a + t + c + g

    Aa = sequence.count('~')
    Tt = sequence.count('!')
    Cc = sequence.count('@')
    Gg = sequence.count('#')
    Duble_masked_count = Aa + Tt + Cc + Gg



    one =  sequence.count('1')
    two =  sequence.count('2')
    three =sequence.count('3')
    four = sequence.count('4')
    none_masked_count = one + two + three + four

    #nusequenced
    N = sequence.count('N')
    n = sequence.count('n')

    #porthin in whole sequence

    unsequenced_region = round((n+N)/Len_seq*100,2)
    if unsequenced_region > 20:
        return 'unsequenced_region > 20'
    one_masked_P = round(((A + T + C + G )/Len_seq*100),2)
    one_masked_C = round(((a+ t + c + g )/Len_seq*100),2)
    Duble_masked = round(((Aa+ Tt + Cc + Gg )/Len_seq*100),2)
    none_masked = round(((one + two + three + four)/Len_seq*100),2)
    masked_portion = one_masked_count_P + one_masked_count_C +Duble_masked_count

    #portions
    try:
        Ap = round(A/one_masked_count_P*100,3)
        Tp = round(T/one_masked_count_P*100,3)
        Cp = round(C/one_masked_count_P*100,3)
        Gp = round(G/one_masked_count_P*100,3)
    except ZeroDivisionError:
        #print('except ZeroDivisionError:','710')
        Ap = 0
        Tp = 0
        Cp = 0
        Gp = 0



    try:
        ap = round(a/one_masked_count_C*100,3)
        tp = round(t/one_masked_count_C*100,3)
        cp = round(c/one_masked_count_C*100,3)
        gp = round(g/one_masked_count_C*100,3)
    except ZeroDivisionError:
        #print('except ZeroDivisionError:',724)
        ap = 0
        tp = 0
        cp = 0
        gp = 0

    try :
        Aap = round(Aa/Duble_masked_count*100,2)
    except ZeroDivisionError:
        #print('except ZeroDivisionError:',733)
        Aap = 0
    try :
        Ttp = round(Tt/Duble_masked_count*100,2)
    except ZeroDivisionError:
        Ttp = 0
    try:
        Ccp = round(Cc/Duble_masked_count*100,2)
    except ZeroDivisionError:
        Ccp = 0
    try:
        Ggp = round(Gg/Duble_masked_count*100,2)
    except ZeroDivisionError:
        Ggp = 0


    try:
        onep = round(one/none_masked_count*100,3)
        twop = round(two/none_masked_count*100,3)
        threep = round(three/none_masked_count*100,3)
        fourp = round(four/none_masked_count*100,3)
    except ZeroDivisionError:
        #print('except ZeroDivisionError:',754)
        onep = 0
        twop = 0
        threep=0
        fourp =0

    #result
    #print('len sequence :',Len_seq)
    #print('one_masked_p',one_masked_P,'one_masked_c',one_masked_C,'%','Duble_masked',Duble_masked,'%','unmasked',none_masked,'%','unsequenced',unsequenced_region,'%')
    #print('Sum of all',one_masked_P + one_masked_C + Duble_masked + none_masked + unsequenced_region)
    #print('all masks portion',round(masked_portion/Len_seq*100,2),'%')


    '''print(r'masks composition')
    print('A',Ap,'%','a',ap,'%','Aa',Aap,'1',onep,'balans result:',round((A+a+Aa)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2))
    print('T',Tp,'%','t',tp,'%','Tt',Ttp,'2',twop,'balans result:',round((T+t+Tt)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2))
    print('C',Cp,'%','c',cp,'%','Cc',Ccp,'3',threep,'balans result:',round((C+c+Cc)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2))
    print('G',Gp,'%','g',gp,'%','Gg',Ggp,'4',fourp,'balans result:',round((G+g+Gg)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2))
    print('sumP:',Ap+Tp+Cp+Gp,'sumC:',ap+tp+cp+gp,'sumPC:',Aap+Ttp+Ccp+Ggp,'sum_unmasked',onep+twop+threep+fourp)'''
    try:
        A_T = round(A/T,3)
        G_C = round(G/C,3)
    except ZeroDivisionError:
        A_T = 0
        G_C = 0

    try:
        a_t =  round(a/t,3)
        g_c =  round(g/c,3)
    except ZeroDivisionError:
        a_t =  0
        g_c =  0

    try:
        Aa_Tt = round(Aa/Tt,3)
        Gg_Cc = round(Gg/Cc,3)
    except:
        Aa_Tt = 0
        Gg_Cc = 0


    try:
        one_two  = round((one/two),3)
        four_three=round((four/three),3)
    except:
        one_two  = 0
        four_three=0
    '''if int(one_masked_C) > 0:
        try:

            print('one_masked_p',one_masked_P,'one_masked_c',one_masked_C,'%','Duble_masked',Duble_masked,'%','unmasked',none_masked,'%','unsequenced',unsequenced_region,'%')

            print('A/T',A_T,'a/t',a_t,'Aa/Tt',Aa_Tt,'A~/T!',one_two , 'GC_countent% :', round((C+G)/(A+T+C+G),2))
            print('G/C',G_C,'g/c',g_c,'Gg/Cc',Gg_Cc,'G#/C@',four_three,'GC_countent% :',round((C+G)/(A+T+C+G),2) )
        except:
            if ZeroDivisionError:
                print('zero error')'''




    #print(r'nucleotide_composition_in_seqenence Pass!!!!')
    # Sample data
    labels = ['+ masked','- masked', '+- masked', 'unmasked_region', 'unsequenced_region']
    sizes = [one_masked_P,one_masked_C,Duble_masked,none_masked,unsequenced_region]# percentages, should add up to 100%
    #print('plot informatio',sizes)
    explode = (0.1, 0.1, 0.1, 0.1,0.1)





    #Plotting the pie chart
    #plt.pie(sizes,labels = labels, explode=explode,  autopct='%1.1f%%', shadow=True, startangle=140)

    # Equal aspect ratio ensures that pie is drawn as a circle
    #plt.axis('equal')

    # Adding a title
    #plt.title('Example Pie Chart')

    # Display the plot


    # N usage for information
    try:
        A_usage = round((A + t + Aa)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2)
        T_usage = round((T + a + Tt)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2)
        C_usage = round((C + g + Cc)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2)
        G_usage = round((G + c + Gg)/(one_masked_count_P+one_masked_count_C+Duble_masked_count)*100,2)
        #print('final nucleotide usage','A',A_usage,'T',T_usage,'C',C_usage,'G',G_usage)
    except ZeroDivisionError:
        #print('except ZeroDivisionError',810)
        A_usage = 0
        T_usage = 0
        C_usage = 0
        G_usage = 0
    #plt.show()
    #print(seq_Nucleotide_composition(sequence) )
    #print(Len_seq)
    #print(sequence)
    #print(A_usage,T_usage,C_usage,G_usage,(round(masked_portion/Len_seq*100,2)))
    return [A_usage,T_usage,C_usage,G_usage,(round(masked_portion/Len_seq*100,2))]