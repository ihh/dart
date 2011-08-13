import sys
from extra_functions import get_cmd_args

usage = '''
Mixture-gen, a script to produce C++ code for mixture branch transducers.  Protpal-specific. 
Usage: mixture_gen.py <number of mixtures>
'''



class Transducer(object):
    def __init__(self, name="Transducer"):
        self.name = name
        self.state_names = []
        self.transitions = {}
        self.transition_weights = {}

    def add_transition(self, state1, state2, weight=1.0):
        for state in [state1, state2]:
            if not state in self.state_names:
                self.state_names.append(state)
        if self.transitions.has_key(state1):
            self.transitions[state1].append(state2)
        else:
            self.transitions[state1] = [state2]
    
        self.transition_weights[(state1, state2)] = weight

    def show_dot(self, weights):
        dot = "digraph %s{\n"%self.name
        for source_state in self.transitions:
            for dest_state in self.transitions[source_state]:
                dot += source_state + "->" + dest_state; 
                if weights:
                    dot += '[label="%s"]'%self.transition_weights[(source_state, dest_state)]
                dot += ';\n'
        print dot+'}'

    def show(self):
        code = '''
// %s - mixture-of-affine gaps branch transducer                                                                      
BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix,
                         double ins_open_rate, double del_open_rate, '''%(num_mixes)
        for i in range(num_mixes):
            code += "double gap_extend_%s, double mix_prior_%s, "%(i,i)
        code = code[:-2] +")\n"

        code +=  \
'''
{
branch_length = branch_length_in;
name = "%s";

vector<sstring> toks = alphabet_in.tokens();
for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
   alphabet.push_back(string(a->c_str()));

alphabet_size = alphabet.size();

double ins_open = 1-exp(-ins_open_rate*branch_length);
double del_open = 1-exp(-del_open_rate*branch_length);

conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length);
'''%self.name


        code += "\n// State identifiers\n"
        for i,state in enumerate(states):
            code += 'int %s = %s;\n'%(state, i)
            
        code += "\n// State list\n"
        for i,state in enumerate(states):
            code += 'states.push_back('+str(i)+');\n'

        code += "\n// State names\n"
        for state in states:
            code += 'state_names.push_back("%s");\n'%state

        code += "\n// State types\n"
        for state in states:
            stype = state_type(state)
            code += 'state_types.push_back("%s");\n'%stype
            
        clearOut = "out.clear();\n"
        clearPair = "transitionPair.clear();\n"
        code += "\n// State transitions\n"
        code += "vector<state> out;\n"
        for source_state in self.transitions:
            code += '\n// Transitions out of ' + source_state + '\n'
            code += clearOut
            for dest_state in self.transitions[source_state]:
                code += "out.push_back(%s);\n"% dest_state
            code += "outgoing[%s] = out;\n"%source_state

            
        code += "\n//Transition weights\n"
        code += "vector<state> transitionPair;\n"
        for source_state in self.transitions:
            code += '\n// Transition weights out of ' + source_state + '\n'
            for dest_state in self.transitions[source_state]:
                code += clearPair
                code += "transitionPair.push_back(%s); transitionPair.push_back(%s);\n"%(source_state, dest_state)
                code += "transition_weight[transitionPair] = %s;\n"%(self.transition_weights[(source_state, dest_state)])
        code += \
'''
vector<double> equilibrium = rate_matrix.create_prior();
for (int i=0; i<states.size(); i++)
    if (state_types[i] == "I") 
        emission_weight_matrix[i] = equilibrium;\n}'''




        print code



def state_type(state):
    if state in ['start','end', 'match']:
        return state[0].upper()
    elif 'wait' in state:
        return 'W'
    elif 'delete' in state:
        return "D"
    elif 'insert' in state:
        return 'I'
    else:
        sys.stderr.write("Unknown state: " + state+'\n')
        sys.exit(1)

num_mixes = int(sys.argv[1])
states = ['start','end','wait_end', 'match', 'wait_match']
for i in range(num_mixes):
    states.append('delete_'+str(i))
    states.append('wait_delete_'+str(i))
    states.append('insert_'+str(i))

branch = Transducer("mixture")

# Start state transitions
for i in range(num_mixes):
    branch.add_transition("start", "insert_%s"%i, "ins_open*mix_prior_%s"%i)
    branch.add_transition("start", "wait_delete_%s"%i, "(1-ins_open)*del_open*mix_prior_%s"%i)
branch.add_transition("start", "wait_match", "(1-ins_open)*(1-del_open)")
branch.add_transition("start", "wait_end", "1.0")


# Insert states' transitions
for i in range(num_mixes):
    branch.add_transition("insert_%s"%i, "insert_%s"%i, "gap_extend_%s"%i)
    for j in range(num_mixes):
        branch.add_transition("insert_%s"%i, "wait_delete_%s"%j, "(1-gap_extend_%s)*del_open*mix_prior_%s"%(i,j))
    branch.add_transition("insert_%s"%i, "wait_match", "(1-gap_extend_%s)*(1-del_open)"%i)
    branch.add_transition("insert_%s"%i, "wait_end", "1.0")
    

# Delete states' transitions
for i in range(num_mixes):
    for j in range(num_mixes):
        branch.add_transition("delete_%s"%i, "insert_%s"%j, "(1-gap_extend_%s)*ins_open*mix_prior_%s"%(i,j))
    branch.add_transition("delete_%s"%i, "wait_match", "(1-gap_extend_%s)*(1-ins_open)"%i)
    branch.add_transition("delete_%s"%i, "wait_delete_%s"%i, "gap_extend_%s"%i)
    branch.add_transition("delete_%s"%i, "wait_end", "1.0")

# Match states' transitions
for i in range(num_mixes):
    branch.add_transition("match", "wait_delete_%s"%i, "del_open*(1-ins_open)*mix_prior_%s"%i)
    branch.add_transition("match", "insert_%s"%i, "ins_open*mix_prior_%s"%i)
branch.add_transition("match", "wait_match", "(1-ins_open)*(1-del_open)")
branch.add_transition("match", "wait_end", "1.0")

# wait_delete  states' transitions
for i in range(num_mixes):
    branch.add_transition("wait_delete_%s"%i, "delete_%s"%i, "1.0")

# wait_match
branch.add_transition("wait_match", "match", "1.0")

# wait_end
branch.add_transition("wait_end", "end", "1.0")
    
if 'cpp' in sys.argv[1:]:
    branch.show()
elif 'dot' in sys.argv[1:]:
    branch.show_dot(weights= ('weights' in sys.argv[1:]))






