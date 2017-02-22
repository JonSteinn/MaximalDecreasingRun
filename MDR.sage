###################################
##### COUNTING AND GENERATING #####
###################################

# input: permutation p
# output: Maximal decreasing run
def MDR(p):
    nextMDR = len(p)
    for i in p:
        if (i == nextMDR):
            nextMDR -= 1
    return len(p) - nextMDR

# input: length of permutations n
# output: [M(n,1),M(n,2),...,M(n,n)]
def distribution_MDR(n):
    M = [0 for i in [1..n]]
    for p in Permutations(n):
        M[MDR(p)-1] += 1
    return M

# input: length of permutations n, MDR k, type of data structure
# output: a list/dictionary with all permutations in S_(n,k)
def generate_MDR(n,data_structure="list"):
    if (data_structure == "list"):
        A = [[] for i in [1..n]]
        for p in Permutations(n):
            A[MDR(p)-1] += [p]
    if (data_structure == "dictionary"):
        A = [{} for i in [1..n]]
        for p in Permutations(n):
            A[MDR(p)-1][p] = 0
    return A

# input: length of permutations n, pattern
# output: [M_{pattern}(n,1),M_{pattern}(n,2),...,M_{pattern}(n,n)]
def distribution_avoiding_MDR(n,pattern):
    M = [0 for i in [1..n]]
    for p in Permutations(n,avoiding=pattern):
        M[MDR(p)-1] += 1
    return M

# input: length of permutation n, MDR k, pattern, type of data structure
# output: a list/dictionary with all permutations in S_{n,k}(pattern)
def generate_avoiding_MDR(n,pattern,data_structure="list"):
    if (data_structure == "list"):
        A = [[] for i in [1..n]]
        for p in Permutations(n,avoiding=pattern):
            A[MDR(p)-1] += [p]
    if (data_structure == "dictionary"):
        A = [{} for i in [1..n]]
        for p in Permutations(n,avoiding=pattern):
            A[MDR(p)-1][p] = 0
    return A

# input: Dyck path d
# output: number of returns to the x-axis
def returns(d):
    c = 0
    height = 0
    for i in d: # d is a list of binary numbers
        if (i):
            height += 1
        else:
            height -= 1
            if (height == 0):
                c += 1
    return c

# input: length of Dyck path n
# output: [|D_{n,1}^{ret}|,|D_{n,2}^{ret}|,...|D_{n,n}^{ret}|]
def distribution_return_Dyck(n):
    M = [0 for i in [1..n]]
    for d in DyckWords(n):
        M[returns(d)-1] += 1
    return M

# input: length of Dyck path n, number of returns k, type of data structure
# output: a list/dictionary of all Dyck paths in D_{n,k}^{ret}
def generate_return_Dyck(n,data_structure="list"):
    if (data_structure == "list"):
        A = [[] for i in [1..n]]
        for p in DyckWords(n):
            A[returns(p)-1] += [p]
    if (data_structure == "dictionary"):
        A = [{} for i in [1..n]]
        for p in DyckWords(n):
            A[returns(p)-1][p] = 0
    return A

# input: Dyck path d
# output: legnth of first descent
def first_descent(d):
    c = 0
    index = 0
    while (d[index] == 1):
        index += 1
        if (index * 2 == len(d)): # first half of length are up steps
            return index
    while (d[index] == 0):
        c += 1
        index += 1
    return c

# input: length of Dyck path n
# output: [|D_{n,1}^{des}|,|D_{n,2}^{des}|,...|D_{n,n}^{des}|]
def distribution_descent_Dyck(n):
    M = [0 for i in [1..n]]
    for d in DyckWords(n):
        M[first_descent(d)-1] += 1
    return M

# input: length of Dyck path n, number of returns k, type of data structure
# output: a list of all Dyck paths in D_{n,k}^{des}
def generate_descent_Dyck(n,data_structure="list"):
    if (data_structure == "list"):
        A = [[] for i in [1..n]]
        for p in DyckWords(n):
            A[first_descent(p)-1] += [p]
    if (data_structure == "dictionary"):
        A = [{} for i in [1..n]]
        for p in DyckWords(n):
            A[first_descent(p)-1][p] = 0
    return A

###################################################
####### Examples: COUNTING AND GENERATING #########
###################################################

#print "Maximal decreasing runs:"
#MDR([1,2,3,4])
#MDR([5,4,3,2,1])
#MDR([6,1,5,2,4,3])
#print

#print "Distribution for permutations of length n with MDR = k:"
#for i in [1..8]:
#    distribution_MDR(i)
#print

#print "Generation of S_{n,k}:"
#generate_MDR(3)
#generate_MDR(3,"dictionary")
#print

#print "Distribution for permutations of length n with MDR = k that avoid a pattern:"
#for i in [1..8]:
#    distribution_avoiding_MDR(i,[1,2,3])
#for i in [1..8]:
#    distribution_avoiding_MDR(i,[3,1,2])
#print

#print "Generation of S_{n,k}(pattern):"
#generate_avoiding_MDR(4,[2,3,1])
#print

#print "Number of returns in a Dyck paths:"
#returns([1,1,1,0,0,0])
#returns([1,0,1,0,1,0])
#returns([1,1,0,1,0,0])
#print

#print "Distribution of Dyck paths with length n and k returns:"
#for i in [1..4]:
#    distribution_return_Dyck(i)
#print

#print "Generation of Dyck paths with length n and k returns:"
#generate_return_Dyck(5)
#print

#print "Length of first descent in a Dyck path:"
#first_descent([1,1,1,1,0,0,0,0])
#first_descent([1,0,1,1,0,1,0,0])
#first_descent([1,1,1,0,1,0,0,0])
#print

#print "Distribution of Dyck paths with length n and first descent of length k:"
#for i in [1..4]:
#    distribution_return_Dyck(i)
#print

#print "Generation of Dyck paths with length n and first descent of length k:"
#generate_descent_Dyck(7,"dictionary")
#print

###############################################
############### recursions gamma ##############
###############################################

# input: a permutation p with length > MDR - 1
# output: a list of permutations with length incremented by 1 by adding 1s
def gamma(p):
    p = [i + 1 for i in p]
    P = [[1] + [i for i in p]] # left end special case
    for i in [0..len(p)-1]: # add 1 to all possible position
        P += [[p[j] for j in [0..i]] + [1] + [p[j] for j in [i+1..len(p)-1]]]
    return P

#################################
####### Examples: gamma #########
#################################

#print "The element corresponding to a permutation in S_{n-1,k} in the image of gamma (must have n < k-1):"
#gamma([4,1,2,3])
#gamma([7,3,2,5,6,1,4])
#print

#################################################
############### recursions r1^(-1) ##############
#################################################

# input: a permutation from S_{n,k}(231)
# output: a permutation p from S_{n-1,k-1}(231) or S_{n,k+1}
def r1_inv(p,k):
    n = len(p)
    if (p[0] == n):
        return Permutation([p[i] for i in [1..n-1]])
    else:
        return Permutation([p[i] for i in [1..p[0]-1]] + [n] + [p[i]-1 for i in [p[0]..n-1]])

##############################
####### Examples: r1 #########
##############################

#print "The element in S_{n-1,k-1}(231) corresponding to a permutation in S_{n,k}(231):"
#r1_inv([1,5,3,2,4],2)
#r1_inv([5,1,3,2,4],2)
#print

#################################################
############### recursions r2^(-1) ##############
#################################################

# input: list of sub permutations L [[...],[...],...,[...]]
# output: same list with elements rotated amongst non empty inner lists
def rotator(L):
    elements = [] # list with elements, no inner lists
    for i in [0..len(L)-1]:
        if (len(L[i]) > 0):
            for j in L[i]:
                elements += [j]
    rotation = {}
    for i in [0..len(elements)-1]: # map each element to it's left neighbour (circular)
        rotation[elements[i]] = elements[(i-1+len(elements)) % len(elements)]
    for i in L: # replace domain elements with corresponding image elements
        if (len(i) > 0):
            for j in [0..len(i)-1]:
                i[j] = rotation.get(i[j])
    return L

# input: 123 avoiding permutation p, MDR k, boolean return_path (for phi_3)
# output: see cases in comments, but always a 123 avoiding permutation
def r2_inv(p,k,return_path=false): # inverse chosen for phi_3, direction dosent matter
    n = len(p)
    if (p[0] == n): # returns a permutation with MDR & length decremented by 1
        if (return_path):
            return (0,Permutation([p[i] for i in [1..n-1]]))
        else:
            return Permutation([p[i] for i in [1..n-1]])
    else: # returns a permutation of same length with MDR incremented by 1
        cut_point = -1 # n-k-1
        last_run_element = -1 # n-k+1
        for i in [0..n-1]: # find n-k-1 and n-k+1
            if p[n-1-i] == n-k+1:
                last_run_element = n-1-i
            if p[n-1-i] == n-k-1:
                cut_point = n-1-i
            if (cut_point >= 0 and last_run_element >= 0):
                break
        if (cut_point < last_run_element): # more trivial subcase, no rotation
            return_value = [p[i] for i in [1..last_run_element]]
            return_value += [n-k] + [p[i] for i in [last_run_element+1..n-1]]
            if (return_path):
                return (1, Permutation(return_value))
            else:
                return Permutation(return_value)
        else:
            # A_{k+1} split in 2 parts (rotated, not rotated)
            B1 = [p[i] for i in [last_run_element+1..cut_point]]
            B2 = [p[i] for i in [cut_point+1..n-1]]
            rotated = []
            temp = []
            next_parse_sign = n
            # gather sub-permutations A_1,A_2,...,A_k to a list
            for i in [1..last_run_element]:
                if (p[i] == next_parse_sign):
                    next_parse_sign -= 1
                    rotated += [temp]
                    temp = []
                else:
                    temp += [p[i]]
            rotated += [B1] # [[A_1],[A_2],...,[A_k], [part of A_{k+1}]]
            rotated = rotator(rotated)
            return_value = []
            run_value = n
             # alternate rotated sub-permutations and run elements
            for i in [0..len(rotated)-1]:
                return_value += rotated[i] + [run_value-i]
            return_value += B2
            if (return_path):
                return (1, Permutation(return_value))
            else:
                return Permutation(return_value)

##############################
####### Examples: r2 #########
##############################

#print "Rotates the non-empty lists within a list by 1 to the right:"
#rotator([[1,2,3],[],[7,4,9]])
#rotator([[1,2,3],[],[7,4,9],[],[],[5],[6]])
#rotator([[1],[2],[3],[4],[5,6]])
#print

#print "The element in S_{n-1,k-1}(123) or S_{n,k+1}(123) corresponding to a permutation in S_{n,k}(231):"
#r2_inv([7,4,6,1,5,3,2],3)
#r2_inv([4,7,1,6,5,3,2],3,true)
#print

##################################################
################### RECURSION r ##################
##################################################

# input: a 132 avoiding permutation p
# output: a 132 avoiding permutation with length & MDR incremented by 1
def r(p,k):
    return Permutation([len(p)+1] + [i for i in p]) # n added to left end

#############################
####### Examples: r #########
#############################

#print "The element in S_{n-1,k-1}(132) corresponding to a permutation in S_{n,k}(132):"
#r([8,7,3,4,5,1,2,6],3)
#r([8,7,3,2,4,5,1,6],3)
#print

##############################
############ PHI 1 ###########
##############################

# input: 132 avoiding permutation p, index i
# output: number of larger elements to the right of i
def larger_to_right(p,i):
    c = 0
    for j in [i+1..len(p)-1]: # iterate right of index
        if (p[j] > p[i]):
            c += 1
        if (c == len(p) - p[i]): # max possible value of c
            return c
    return c

# input: 132 avoiding permutation p
# output: Dyck path
def Krattenthaler(p):
    """
    Krattenthaler's algorithm:
        - input: permutation p, outout: Dyck path d
        - Read through p from left to right. For element i, count number of
          larger elements to it's right. Suppose there are x elements larger
          than i to it's right. If x is larger or equal than current height
          of d (as it stands now) then add (to right end of d) a as many up
          steps as needed to d so that we can take a last step (for i) from
          height x+1 to x. If x is less than current height than as many
          down steps are added as needed (can be none) to achieve the same.
    """
    d = []
    curr_height = 0
    for i in [0..len(p)-1]:
        r = larger_to_right(p,i)
        while (curr_height <= r):
            curr_height += 1
            d += [1]
        while (curr_height > r):
            curr_height -= 1
            d += [0]
    return d

# input: 231 avoiding permutation p
# output: Dyck path with p's MDR returns
def phi_1(p,k):
    d = [] # stores return value (Dyck path)
    sub = [] # stores sub permutation for helper maps
    nextMDR = len(p) # next run element
    for i in [0..len(p)-1]:
        if (p[i] != nextMDR): # add non-run elements to sub permutation
            sub += [p[i]]
        else: # (p[i] == nextMDR)
            nextMDR -= 1 # update next run element
            if (len(sub) == 0): # empty sub permutation case
                d += [1,0]
            else:
                # reverse: 231->132, also handles lovering elements to 1,2,...,x
                temp = Permutation([j - (min(sub)-1) for j in sub]).reverse()
                d += ([1] + Krattenthaler(temp) + [0])
            sub = []
    return DyckWord(d)

#################################
####### Examples: phi_1 #########
#################################

#print "Number of elements to the right of index i that are larger than the element in i:"
#larger_to_right([1,2,3,4,5,6], 0)
#larger_to_right([6,5,4,3,2,1], 0)
#larger_to_right([1,7,5,4,3,8,6,2], 3)
#print


#print "Krattenthaler's bijection from 132 avoiding permutation to Dyck paths:"
#Krattenthaler([6,4,3,1,2,5,7,8])
#Krattenthaler([8,5,6,7,3,1,2,4])
#print

#print "The element in D_{n,k}^{ret} corresponding to a permutation in S_{n,k}(231):"
#phi_1([5,3,1,2,4,9,8,7,6],4)
#phi_1([1,4,3,2,5,9,8,7,6],4)
#print

#####################################
############## PHI 2 ################
#####################################

# input: 213 avoiding permutation p, MDR k
# output: [A_1,A_2,...,A_(k+1)] (see structure of 213's)
def parser(p,k):
    n = len(p)
    parsed = [] # hols all sub permutations, return value
    temp = [] # holds each sub permutation
    nextMDR = n # next run value
    for i in [0..n-1]:
        if (p[i] == nextMDR): # if in {n,n-1,...,n-k+1}
            parsed += [temp]
            temp = []
            if (nextMDR == n-k+1): # end case
                parsed += [[p[j] for j in [i+1..n-1]]]
                break
            nextMDR -= 1
        else:
            temp += [p[i]]
    return parsed

# input [A_1,A_2,...,A_(k+1)] (parsed)
# output: [alpha_0,alpha_1,...,alphak] (max values of intervals, see text)
def alphas(parsed):
    alph = [0] # alpha_0 = 0 for all parsings
    for i in [1..len(parsed)-1]:
        if (len(parsed[i-1]) == 0): # if empty sub permutation
            alph += [alph[i-1]] # alpha_i copies alpha_{i-1}
        else:
            alph += [max(parsed[i-1])] # largests element of sub permutation
    return alph

# input: sub permutation A_(k+1), alpha-intervals
# output: [A^1, A^2,...,A^k] (reverse order as in A_{k+1})
def parser_k_1(sub, alph):
    parts = [[] for i in [1..len(alph)-1]] # parsin of A_{k+1} into A^k...A^1
    temp = [] # stores sub permutation of a sub permutation
    interval = 0 # which interval
    while (interval < len(alph) - 1):
        if (alph[interval + 1] != alph[interval]): # if interval is not empty
            for i in sub: # iterate through A_{k+1}, store all values on interval
                if (alph[interval] < i and alph[interval +  1] > i):
                    temp += [i]
            parts[interval] += temp # add to list of lists
        interval += 1
        temp = []
    return parts

# input: 213 avoiding permutations p, MDR k
# output: g of the permutation (see text)
def g(p,k):
    parsedLeft = parser(p, k) # A_1, A_2, ..., A_k (and A_{k+1})
    # A^1, A^2, ..., A_k
    parsedRight = parser_k_1(parsedLeft[len(parsedLeft)-1], alphas(parsedLeft))
    # [[A_1A^1],[A_2A^2],...,[A_kA^k]]
    permutations = [parsedLeft[i] + parsedRight[i] for i in [0..k-1]]
    return permutations

# input: 213 avoiding permutation p
# output: 231 avoiding permutation
def omega2(p):
    if (len(p) == 0):
        return []
    # decrement to a permutation on {1,2,...,|p|}
    m = min(p)
    p = Permutation([i - (m-1) for i in p])
    p = p.reverse().complement().reverse() # composition of symmetry maps
    p = [i + (m-1) for i in p] # increment back to original values
    return p

# input: 213 avoiding permutation p, MDR k
# output: 231 avoiding permutation with k MDR
def phi_2(p,k):
    n = len(p)
    parts = g(p,k)
    return_permutation = []
    for i in [0..len(parts)-1]:
        return_permutation += omega2(parts[i]) + [n-i]
    return Permutation(return_permutation)

#################################
####### Examples: phi_2 #########
#################################

#print "Parses a 213 avoiding permutation with k MDR into A_1, A_2,...,A_{k+1}:"
#parser([9,8,7,6,4,5,1,3,2],5)
#parser([5,6,8,7,4,1,2,3],2)
#print

#print "Alpha intervals of a [A_1,A_2,...,A_k] parsing of a 213 avoiding permutation:"
#alphas(parser([9,8,7,6,4,5,1,3,2],5))
#alphas([[2],[],[],[],[4],[3,1]])
#print

#print "Parses A_{k+1} (of the previous parsing) according to alpha intervals:"
#parser_k_1([3,1],[[2],[],[],[],[],[4],[3,1]])
#parser_k_1(parser([8,3,6,7,4,5,2,1],MDR([8,3,6,7,4,5,2,1]))[len(parser([8,3,6,7,4,5,2,1],MDR([8,3,6,7,4,5,2,1])))-1],alphas(parser([8,3,6,7,4,5,2,1],MDR([8,3,6,7,4,5,2,1]))))
#print

#print "Corresponding elements in the image of g for a permutation p in S_{n,k}(213):"
#g([2,3,9,8,7,6,4,5,1],5)
#g([1,2,6,8,7,3,5,4],2)
#print

#print "The symmetry map y=rev(com(rev(x))) from 213 avoiding permutations to 231 avoiding:"
#omega2([3,10,9,8,7,5,6,4,2,1])
#omega2([12,16,17,15,13,14,11])
#print

#print "The element of S_{n,k}(231) corresponding to a permutation in S_{n,k}(213)"
#phi_2([6,7,1,3,5,4,2],1)
#phi_2([10,4,9,8,5,7,6,3,2,1],5)
#print

##############################################
################## Phi 3 #####################
##############################################

# input: a permutation p from S_{n-1,k-1}(231)
# output: a permutation from S_{n,k}(231)
def first_case_r1(p):
    return Permutation([len(p) + 1] + [i for i in p]) # add n to left end

# input: a permutation p from S_{n,k+1}(231)
# output: a permutation from S_{n,k}(231)
def second_case_r1(p):
    temp = []
    index = 1
    while (1):
        if (p[index-1] == len(p)): # find index of n
            break
        else:
            temp += [p[index-1]] # elements up to n
            index += 1
    return [index] + temp + [p[i-1] + 1 for i in [index+1..len(p)]]

# input: a 123 avoiding permutation
# output: a 231 avoiding permutation
def phi_3(p,k,return_path=false):
    permutation_path = [p]
    path = [] # stores unique path from (n,k) to (1,1)
    while (p != [1]): # map p (and store in self) until we reach (1,1)
        A = r2_inv(p,MDR(p),true) # A = [path value, image permutation]
        path = [A[0]] + path # add to front to get path from (1,1) to (n,k)
        p = A[1]
        permutation_path += [p]
    for i in path: # follow the same path for r1
        if (i):
            p = second_case_r1(p)
        else:
            p = first_case_r1(p)
        permutation_path += [p]
    if (return_path):
        return [Permutation(p),permutation_path]
    else:
        return Permutation(p)

#################################
####### Examples: phi_3 #########
#################################

#print "The element of S_{n,k}(231) corresponding to a permutation in S_{n,k}(123):"
#phi_3([2,3,1],1,true)
#phi_3([6,4,2,5,1,3],2)
#phi_3([8,4,7,3,6,5,1,2],4)
#print

##########################################################
################### phi 4 ################################
##########################################################

# input: a 132 avoiding permutation, MDR k
# output: a Dyck path with first descent of length MDR
def phi_4(p,k):
    n = len(p)
    if (n == k): # sepcial case, [u,u,...,u,d,d,...,d]
        return DyckWord([1 for i in [1..k]] + [0 for i in [1..k]])
    A = [[],[]] # [A_1,A_2]
    index = 0
    for i in [k-1..n-1]: # iterate through permutation (can skip a few)
        if (p[i] == n-k+1): # n-k+1
            index += 1 # start adding to A[1]
        else:
            A[index] += [p[i]] # non run elements
    if (len(A[0]) > 0): # if A_1 is not empty
        m = min(A[0])
        if (m != 1): # decrement to permutation on {1,2,...,|A_1|}
            A[0] = [i-(m-1) for i in A[0]]
    D = [Krattenthaler(A[1]),Krattenthaler(A[0])] # notice order of A's
    d = []
    index = 0
    if (len(D[0]) > 0):
        while (D[0][index]): # add consecutive up steps
            d += [1]
            index += 1
        d += [1 for i in [0..k-1]] # k up steps
        d += [0 for i in [0..k-1]] # k down steps
        d += [i for i in D[1]]
        d += [D[0][i] for i in [index..len(D[0])-1]] # add remaining steps
        return DyckWord(d)
    else:
        return DyckWord([1 for i in [1..k]] + [0 for i in [1..k]] + [i for i in D[1]])

#################################
####### Examples: phi_4 #########
#################################

#print "The element of D_{n,k}^{des} corresponding to a permutation in S_{n,k}(132):"
#phi_4([10,9,8,7,5,4,6,3,2,1],5)
#phi_4([6,7,5,3,2,4,1],5)
#print

#################################################
###################### phi_5 ####################
#################################################

# input: 312 avoiding permutation p, MDR k
# output: h of that permutation (see text)
def h(p,k):
    n = len(p)
    A2 = []
    a2_max = 0
    passed = false
    for i in [k..n-1]: # iterate through p to find A_2 (can skip some)
        if (p[i] == n - k + 1): # everything after n-k+1 belongs to A_2
            passed = true
            continue
        if (passed):
            A2 += [p[i]]
            if (a2_max < p[i]): # finds max of A_2
                a2_max = p[i]
    A1 = []
    for i in [0..n-1]: # iterate through p to find A_1
        if (p[i] == n): # won't go further than n
            break
        A1 += [p[i]]
    temp = [[],[]] # stores A_1 when split by size
    for i in [0..len(A1)-1]: # iteration A1
        if (A1[i] > a2_max): # larger than max(A2) elements
            temp[1] += [A1[j] for j in [i..len(A1)-1]]
            break
        temp[0] += [A1[i]] # smaller than max(A2) elements
    return [temp[1],temp[0]+A2]

# input: 312 avoiding permutation on {a,...,b}
# output: complement of that permutation (132 avoiding)
def com(p):
    if (len(p) == 0): # special case
        return []
    m = min(p)
    p = [i-(m-1) for i in p] # decrement to permutation on {1,2,...,len(p)}
    p = Permutation(p)
    p = p.complement() # 312 -> 132
    p = [i+(m-1) for i in p] #increment to original
    return p

# input: 312 avoiding permutation
# output: 132 avoiding permutation
def phi_5(p,k):
    temp = h(p,k)
    n = len(p)
    perm = []
    for i in [0..k-2]: # add first k-1 run elements
        perm += [n-i]
    perm += com(temp[0])
    perm += [n-k+1] # add last run element
    perm += com(temp[1])
    return Permutation(perm)

#################################
####### Examples: phi_5 #########
#################################

#print "Corresponding elements in the image of h for a permutation p in S_{n,k}(312):"
#h([3,4,2,1,5],1)
#h([1,3,4,5,8,7,6,2],3)
#print

#print "Symmetry map y = com(x), 312->132 avoiding"
#com([1,3,5,4,2])
#com([8,7,10,9])
#print

#print "Element in S_{n,k}(132) corresponding to a permutation in S_{n,k}(312):"
#phi_5([1,2,3,7,6,5,8,4,10,9],2)
#phi_5([4,5,3,2,12,11,10,9,8,7,6,1],7)
#print

###################################################
############### Validation of maps ################
###################################################

def identity_checker(n):
    for i in [1..n]:
        P = generate_MDR(i)
        P1 = generate_avoiding_MDR(i,[2,3,1])
        P2 = generate_avoiding_MDR(i,[2,1,3])
        P3 = generate_avoiding_MDR(i,[1,2,3])
        P4 = generate_avoiding_MDR(i,[3,1,2])
        P5 = generate_avoiding_MDR(i,[1,3,2])
        D1 = generate_return_Dyck(i)
        D2 = generate_descent_Dyck(i)
        for j in [1..i]:
            A = [len(P[j-1])] # lengths of S_{i,j}
            if (j < i):
                A += [sum([x*(binomial(x-1,j-1))*factorial(i-j-1) for x in [0..i-1]])]
                A += [j * factorial(i) / (factorial(j+1))]
            else:
                A += [1]
            for x in [0..len(A)-2]: # are all values in list equal?
                if (A[x] != A[x+1]):
                    print "identiy does not hold"
            B = [len(P1[j-1]),len(P2[j-1]),len(P3[j-1]),len(D1[j-1])]
            B += [binomial(2*i - j - 1, i - j) * j / i]
            for x in [0..len(B)-2]: # are all values in list equal?
                if (B[x] != B[x+1]):
                    print "identity does not hold"
            C = [len(P4[j-1]),len(P5[j-1]),len(D2[j-1])]
            if (i > j):
                C += [catalan_number(i-j+1)-catalan_number(i-j)]
            else:
                C += [1]
            for x in [0..len(C)-2]: # are all values in list equal?
                if (C[x] != C[x+1]):
                    print "identity does not hold"
            print i, j, "DONE"

# A general purpose map checker!
# Checks if the map func is a bijection from domain to image
# Dommain must be a list of permutations, image must be a hashmap
# MDR k is needed for some function, might be unused for others
def map_checker(domain, image, func, k):
    A = {}
    for p in domain:
        p = func(p,k)
        if (A.has_key(p)): # have we mapped to this element before?
            print "Not injective!"
        A[p] = 0
        if (not image.has_key(p)): # is this element in the image?
            print "Invalid mapping!"
    if (len(A) != len(image)): # is there anything we missed?
        print "Not surjective!"
    return "DONE!"

def check_gamma(n):
    for i in [1..n]:
        X = generate_MDR(i)
        Y = generate_MDR(i+1,"dictionary")
        for j in [1..i-1]:
            A = {}
            dom = []
            for p in X[j-1]:
                dom += gamma(p)
            print i, j, map_checker(dom, Y[j-1],Permutation,j)

def check_r1_inv(n):
    Y = [generate_avoiding_MDR(i,[2,3,1],"dictionary") for i in [1..n]]
    for i in [2..n]:
        X = generate_avoiding_MDR(i,[2,3,1])
        for j in [1..i]:
            A = {}
            IM = {}
            if (j < i):
                IM.update(Y[i-1][j])
            if (j > 1):
                IM.update(Y[i-2][j-2])
            print i,j,map_checker(X[j-1],IM,r1_inv,j)

def check_r2_inv(n):
    Y = [generate_avoiding_MDR(i,[1,2,3],"dictionary") for i in [1..n]]
    for i in [2..n]:
        X = generate_avoiding_MDR(i,[1,2,3])
        for j in [1..i]:
            A = {}
            IM = {}
            if (j < i):
                IM.update(Y[i-1][j])
            if (j > 1):
                IM.update(Y[i-2][j-2])
            print i,j,map_checker(X[j-1],IM,r2_inv,j)

def check_r(n):
    X = [generate_avoiding_MDR(i,[1,3,2]) for i in [1..n]]
    Y = [generate_avoiding_MDR(i,[1,3,2],"dictionary") for i in [1..n+1]]
    for i in [1..n]:
        for j in [1..i]:
            print i,j,map_checker(X[i-1][j-1],Y[i][j],r,j)

def check_phi_1(n):
    for i in [1..n]:
        X = generate_avoiding_MDR(i,[2,3,1])
        Y = generate_return_Dyck(i,"dictionary")
        for j in [1..i]:
            print i,j,map_checker(X[j-1],Y[j-1],phi_1,j)

def check_phi_2(n):
    for i in [1..n]:
        X = generate_avoiding_MDR(i,[2,1,3])
        Y = generate_avoiding_MDR(i,[2,3,1],"dictionary")
        for j in [1..i]:
            print i,j,map_checker(X[j-1],Y[j-1],phi_2,j)

def check_phi_3(n):
    for i in [1..n]:
        X = generate_avoiding_MDR(i,[1,2,3])
        Y = generate_avoiding_MDR(i,[2,3,1],"dictionary")
        for j in [1..i]:
            print i,j,map_checker(X[j-1],Y[j-1],phi_3,j)

def check_phi_4(n):
    for i in [1..n]:
        X = generate_avoiding_MDR(i,[1,3,2])
        Y = generate_descent_Dyck(i,"dictionary")
        for j in [1..i]:
            print i,j,map_checker(X[j-1],Y[j-1],phi_4,j)

def check_phi_5(n):
    for i in [1..n]:
        X = generate_avoiding_MDR(i,[3,1,2])
        Y = generate_avoiding_MDR(i,[1,3,2],"dictionary")
        for j in [1..i]:
            print i,j,map_checker(X[j-1],Y[j-1],phi_5,j)

#######################################
####### Examples: Map checker #########
#######################################

#print "Gamma:"
#%time check_gamma(10)
#print "r1:"
#%time check_r1_inv(10)
#print "r2:"
#%time check_r2_inv(10)
#print "r:"
#%time check_r(10)
#print "phi1:"
#%time check_phi_1(10)
#print "phi2:"
#%time check_phi_2(10)
#print "phi3:"
#%time check_phi_3(10)
#print "phi4:"
#%time check_phi_4(10)
#print "phi5:"
#%time check_phi_5(8)

###################################################
############## Demonstration of maps ##############
###################################################

def show_gamma(n,k):
    for p in generate_MDR(n)[k-1]:
        print p, " -> ", gamma(p)

def show_r1_inv(n,k):
    for p in generate_avoiding_MDR(n,[2,3,1])[k-1]:
        print p, " -> ", r1_inv(p,k)

def show_r2_inv(n,k):
    for p in generate_avoiding_MDR(n,[1,2,3])[k-1]:
        print p, " -> ", r2_inv(p,k)

def show_r(n,k):
    for p in generate_avoiding_MDR(n,[1,3,2])[k-1]:
        print p, " -> ", r(p,k)

# input: dyck path d
# output: list of coordinates it passes through in correct order
def Dyck_to_cords(d):
    current_cord = [0,0]
    cords = [(current_cord[0],current_cord[1])]
    for i in d:
        current_cord[0] += 1
        if (i):
            current_cord[1] += 1
        else:
            current_cord[1] -= 1
        cords += [(current_cord[0],current_cord[1])]
    return cords

# print/plot elements of domain with corresponding elements of image
def show_phi_1(n,k,show_plot=false):
    if (show_plot):
        for p in generate_avoiding_MDR(n,[2,3,1])[k-1]:
            print p
            cords = Dyck_to_cords(phi_1(p,k))
            lines = line2d(cords)
            points = point2d(cords,size=50,color='red')
            plot(lines+points).show(figsize=3)
    else:
        for p in generate_avoiding_MDR(n,[2,3,1])[k-1]:
            print p, " -> ", phi_1(p,k)

def show_phi_2(n,k):
    for p in generate_avoiding_MDR(n,[2,1,3],"list")[k-1]:
        print p, " -> ", phi_2(p,k)

# print elements of domain with corresponding elements of image
def show_phi_3(n,k,show_path=false):
    if (show_path):
        for p in generate_avoiding_MDR(n,[1,2,3])[k-1]:
            print [i for i in phi_3(p,k,true)[1]]
    else:
        for p in generate_avoiding_MDR(n,[1,2,3])[k-1]:
            print p, " -> ", phi_3(p,k)

def show_phi_4(n,k,show_plot=false):
    if (show_plot):
        for p in generate_avoiding_MDR(n,[1,3,2])[k-1]:
            print p
            cords = Dyck_to_cords(phi_4(p,k))
            lines = line2d(cords)
            points = point2d(cords,size=50,color='red')
            plot(lines+points).show(figsize=3)
    else:
        for p in generate_avoiding_MDR(n,[1,3,2])[k-1]:
            print p, " -> ", phi_4(p,k)

def show_phi_5(n,k):
    for p in generate_avoiding_MDR(n,[3,1,2])[k-1]:
        print p, " -> ", phi_5(p,k)

#############################################
####### Examples: Map demonstration #########
#############################################

#print "gamma:"
#show_gamma(5,4)
#print "r1:"
#show_r1_inv(7,4)
#print "r2:"
#show_r2_inv(5,3)
#print "r:"
#show_r(10,6)
#Dyck_to_cords([1,1,0,0,1,1,0,1,0,0])
#print "phi1"
#show_phi_1(3,1)
#show_phi_1(9,7,true)
#print "phi2"
#show_phi_2(5,3)
#print "phi3"
#show_phi_3(5,1)
#show_phi_3(4,2,true)
#print "phi4"
#show_phi_4(4,2)
#show_phi_4(11,8,true)
#print "phi5"
#show_phi_5(7,3)
