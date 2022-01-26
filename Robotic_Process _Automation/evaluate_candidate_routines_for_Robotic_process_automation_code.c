
/* Program to evaluate candidate routines for Robotic Process Automation.

 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>

#define MAXCHARS 2000    /* max chars per word */
#define INITIAL 100    /* initial size of word array */
#define ASIZE 26    
#define UNSIGNED 256
#define ON 1
#define OFF 0
#define TURNOFF 2
/* type definitions ----------------------------------------------------------*/

 // state (values of the 26 Boolean variables)
 typedef unsigned char state_t[UNSIGNED];

 // action
 typedef struct action action_t;
 struct action {
     char name;        // action name
     state_t precon;   // precondition
     state_t effect;   // effect
 };

 // step in a trace
 typedef struct step step_t;
 struct step {
     action_t *action; // pointer to an action performed at this step
     step_t   *next;   // pointer to the next step in this trace
 };

 // trace (implemented as a linked list)
 typedef struct {
     step_t *head;     // pointer to the step in the head of the trace
     step_t *tail;     // pointer to the step in the tail of the trace
 } trace_t;

typedef char word_t[MAXCHARS+1];

int getword(word_t W, int limit);
void exit_if_null(void *ptr, char *msg);
trace_t* make_empty_trace(void);
trace_t* insert_at_tail(trace_t*, action_t*);
void free_trace(trace_t*);
void print_int_array(char str[], int start, int end);
void build_state(char* str, state_t state);
void turn_off(char* str, state_t* state);
action_t* read_action(char str[], int start, int end);
int check_precon(state_t pre, state_t state);
void do_effect(state_t *effect, state_t *state);
void exit_if_null(void *ptr, char *msg);
int is_same(state_t A, state_t B);
void copy_states(state_t source, state_t dest);
int is_modified(state_t effect, state_t state);
void print_out_state(state_t state);
int search_for_stageone(trace_t* trace, char* commands, int cm, 
                        state_t variables, int start_pos);
int search_for_stagetwo(trace_t* trace, state_t copy_state, char* commands, 
                        int cm, int start_pos, state_t variables_state);
int find_distinct_actions(char* actions);
void update_variables(state_t* state, state_t* effect);
void  turn_it_off(state_t *effect, state_t *state);


int
main(int argc, char *argv[]) {
	word_t one_word;
	char **all_words;
	size_t current_size=INITIAL;
	int numdistinct=0, i, j;
    state_t initial_state, new_state, test_state, second_copy, variables_state;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    action_t** all_actions =(action_t**)malloc(INITIAL*sizeof(*all_actions));   
    /*an array of action_t*/
    exit_if_null(all_actions, "initial allocation");
    trace_t* command;
    step_t* tmp;
    int distinct_action=0, input_len, status, m, numvalid=0, cand_start=0, 
    stagetwo_start, actionstart, numaction=0;

    all_words = (char**)malloc(INITIAL*sizeof(*all_words));
	exit_if_null(all_words, "initial allocation");

	while (getword(one_word, MAXCHARS) != EOF) {
        
		/* a new word exists, check if there is space? */
		if (numdistinct == current_size) {
				current_size *= 2;
				all_words = realloc(all_words, current_size*sizeof(*all_words));
				exit_if_null(all_words, "reallocation");
		}
		/* there is definitely space in array */
		all_words[numdistinct] = (char*)malloc(1+strlen(one_word));
		exit_if_null(all_words[numdistinct], "string malloc");
		/* and there is also a space for the new word */
		strcpy(all_words[numdistinct], one_word);
		numdistinct += 1;
	}

    // initialize the state and build the initial state
    for (i=0;i<ASIZE;i++) {
        initial_state[(int)alpha[i]] = OFF;
    }
    for (i=0;all_words[0][i]!='\0';i++){
        if (isalpha(all_words[0][i])) {
            initial_state[(int)all_words[0][i]] = ON;
        }
    }
   


    actionstart = 0;
    // build the input trace, store each action into an array 
    while (actionstart < strlen(all_words[1])) {
        for (i=actionstart;all_words[1][i]!='\0';i++) {
            if (all_words[1][i] == '\n') {
                if (numaction *sizeof(action_t*) == current_size) {
			        current_size *= 2;
				    all_actions = realloc(all_actions, 
                    current_size*sizeof(*all_actions));
			 	    exit_if_null(all_actions, "reallocation");
			        }
                all_actions[numaction] =(action_t*)malloc(sizeof(action_t));
		        exit_if_null(all_actions[numaction],"string malloc");
                all_actions[numaction] = read_action(all_words[1], actionstart, 
                                                     i-1); 
                numaction++;
                break;
            }     
        }
        actionstart = i+1;
    }


    // build the command trace into a linked list
    command = make_empty_trace();
    for (i=0;all_words[2][i]!='\n';i++) {
        for (j=0;j<numaction;j++) {
            if (all_actions[j]->name == all_words[2][i]) {
                command = insert_at_tail(command, all_actions[j]);
            }
        }
    }

    // find distinct actions 
    if (all_words[2]){
        distinct_action = find_distinct_actions(all_words[2]);
    }


    // copy a state - test_state for testing each action is valid
    copy_states(initial_state, test_state);

    printf("==STAGE 0===============================\n");
    printf("Number of distinct actions: %d\n", distinct_action);
    input_len = i-1;
    printf("Length of the input trace: %d\n", input_len);

    copy_states(initial_state, test_state);
    // check if the trace is valid 
    tmp = command->head;
    
    /*transverse through the list, see if conditions are valid*/
    while (tmp!=NULL) {
        status = check_precon(tmp->action->precon, test_state);
        if (status != 1) {
            break;
        }
        if (status == 1) numvalid++;
        do_effect(&tmp->action->effect, &test_state);
        tmp = tmp->next;  
    }
    if (status == 0) {
       printf("Trace status: invalid\n");
    }
    else {
        printf("Trace status: valid\n");
    }
    
    printf("----------------------------------------\n");
    printf("  %s\n", alpha);
    printf("> ");
    print_out_state(initial_state);


    copy_states(initial_state, new_state);
    tmp = command->head;
    /*transverse through the list, and print out updated states*/
    for (i=0;i<numvalid;i++) {
        printf("%c ", tmp->action->name);
        do_effect(&tmp->action->effect, &new_state);
        print_out_state(new_state);
        tmp = tmp->next;  
    }
 
    
    // if there is no stage 1, or trace invalid, then exit
    if (numdistinct == 3 || (status == 0)) {
        free_trace(command);
        for (i=0; i<numdistinct; i++) {
		    free(all_words[i]);
		    all_words[i] = NULL;
	    }
	    free(all_words);
	    all_words = NULL;
        for (i=0; i<numaction; i++) {
		free(all_actions[i]);
		all_actions[i] = NULL;
	}
	free(all_actions);
	all_actions = NULL;
        return EXIT_SUCCESS;
    }

    printf("==STAGE 1===============================\n");
     // initialize the variables state
    for (i=0;i<ASIZE;i++) {
        variables_state[(int)alpha[i]] = OFF;
    }
    while (cand_start < strlen(all_words[3])) {
    // build an expected array of variables that are turned on and off! 
    // we build an expected effect state 
        for (i=cand_start; all_words[3][i]!='\n';i++) {
            for (j=0;j < numaction;j++) {
                if (all_actions[j]->name == all_words[3][i]) {
                    update_variables(&variables_state, &all_actions[j]->effect);
                }
            }
        }
    
    /*now, search through the trace 'command' and find the same state effect*/
    printf("Candidate routine: ");
    print_int_array(all_words[3], cand_start, i);
        
    // to prevent overlapping, we skip any subsequence that was searched before
        m=0;
        while (m < input_len) {
            m = search_for_stageone(command, all_words[2], input_len, 
            variables_state, m);
        }
        if (all_words[3][i+1]!='\0') {
            printf("----------------------------------------\n");
        }
        for (j=0;j<ASIZE;j++) {
            variables_state[(int)alpha[j]] = OFF;
        }
        cand_start = i+1;
    }

    // if there is no stage 2, exit
    if (numdistinct == 4) {
        free_trace(command);
        for (i=0; i<numdistinct; i++) {
		    free(all_words[i]);
		    all_words[i] = NULL;
	    }
	    free(all_words);
	    all_words = NULL;
        for (i=0; i<numaction; i++) {
		    free(all_actions[i]);
		    all_actions[i] = NULL;
	    }
	    free(all_actions);
	    all_actions = NULL;
        return EXIT_SUCCESS;
    }

    // for stage 2

    /*intialize a second state*/

    for (i=0;i<ASIZE;i++) {
        second_copy[(int)alpha[i]] = OFF;
    }

    for (i=0;i<ASIZE;i++) {
        variables_state[(int)alpha[i]] = OFF;
    }

    printf("==STAGE 2===============================\n");
    stagetwo_start = 0;
    while (stagetwo_start < strlen(all_words[4])) {
    
        // we build an expected effect state
        for(i=stagetwo_start; all_words[4][i]!='\n';i++) {
            for (j=0;j<numaction;j++) {
                if (all_actions[j]->name == all_words[4][i]) {
                    do_effect(&all_actions[j]->effect, &second_copy);
                    update_variables(&variables_state, &all_actions[j]->effect);
                }
            }
        }
  
        /*now, search through the trace 'command' and find the same effect*/
        /*search for effect*/ 
        printf("Candidate routine: ");
        print_int_array(all_words[4], stagetwo_start, i);
        // to prevent overlapping! 
        m=0;
        while (m<input_len) {
            m = search_for_stagetwo(command, second_copy, all_words[2], 
            input_len, m, variables_state);
        }
        if (all_words[4][i+1]!='\0'){
            printf("----------------------------------------\n");
        }

        for (j=0;j<ASIZE;j++) {
            second_copy[(int)alpha[j]] = OFF;
        }
        for (j=0;j<ASIZE;j++) {
            variables_state[(int)alpha[j]] = OFF;
        }
        stagetwo_start = i+1;
    }

    printf("==THE END===============================\n");

    free_trace(command);
    for (i=0; i<numdistinct; i++) {
		free(all_words[i]);
		all_words[i] = NULL;
	}
	free(all_words);
	all_words = NULL;
    for (i=0; i<numaction; i++) {
		free(all_actions[i]);
		all_actions[i] = NULL;
	}
	free(all_actions);
	all_actions = NULL;
	return EXIT_SUCCESS;
}


int
getword(char W[], int limit) {
	int c, len=0;

	/* first, skip over any non alphabetics */
	while ((c=getchar())!=EOF && !isalpha(c)) {
        if (c == '#'){
            W[len] = '\0';
            return 0;
        }
        // in case there is no preconditions for a given command
        if (c == ':'){
            while (len<limit && (c=getchar())!=EOF && c!='#') {
		    /* another character to be stored */
                W[len++] = c;
	    }
	    /* now close off the string */
	    W[len] = '\0';
	    return 0;
        }     
		/* do nothing more */
	}
	if (c==EOF) {
		return EOF;
	}
	/* ok, first character of next word has been found */
	W[len++] = c;
	while (len<limit && (c=getchar())!=EOF && c!='#') {
		/* another character to be stored */
            W[len++] = c;
	}
	/* now close off the string */
	W[len] = '\0';
	return 0;
}

/* =====================================================================
   Program written by Alistair Moffat, 'getword.c'
   as an example in the lecture slide
   modified by Na Zhao

   See http://people.eng.unimelb.edu.au/ammoffat/ppsaa/ for further
   information.
   ================================================================== */

int 
find_distinct_actions(char* actions) {
	int i, distinct=0;
    state_t wordfreq;
    for (i=0;i<UNSIGNED;i++){
        wordfreq[i] = 0;
    }
    for (i=0;actions[i]!='\0';i++){
        if (isalpha(actions[i])) {
            wordfreq[(int)actions[i]] = 1;
        }
    }
    for (i=0;i<UNSIGNED;i++) {
        if (wordfreq[i] == 1){
            distinct++;
        }
    }
	return distinct;
}
/*this function implements stage one, it searches the input trace and find 
the shorted string that 
do the correct effect*/
int 
search_for_stageone(trace_t* trace, char* commands, int cm, state_t variables, 
                    int start_pos) {
    step_t* tmp;
    state_t off;
    state_t new;
    int i=0, k=0, j;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";

    assert(trace->head!=NULL);


    for (j=0;j<ASIZE;j++){
        off[(int)alpha[j]] = OFF;
    }
    for (j=0;j<ASIZE;j++){
        if (variables[(int)alpha[j]] == TURNOFF){
            new[(int)alpha[j]] = 1;
        }
        else {
            new[(int)alpha[j]] = variables[(int)alpha[j]];
        }
    }
     
    tmp = trace->head;
    for (i=0;i<start_pos;i++){
        tmp = tmp->next;
    }
    
    /*step through each command in trace, then do effect, check if same*/
    while (tmp != NULL) { 

        // need to store variables turned on and off into an array, check if 
        // those variables are turned on or off!
     
        if (is_modified(tmp->action->effect, new)) {
            return start_pos+1;
        }
        // then update the effect
            turn_it_off(&tmp->action->effect, &off);
            k++;
            if (is_same(off, variables)) {
                printf("%5d: ", start_pos);
                print_int_array(commands, start_pos, start_pos+k); 
                return start_pos+k;
            }
          
            tmp = tmp->next;
    }   
    return start_pos+1;
}


void update_variables(state_t* state, state_t* effect) {
    int i;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    for (i=0;i<ASIZE;i++) {
        if ((*effect)[(int)alpha[i]]!= OFF) {
            (*state)[(int)alpha[i]] = (*effect)[(int)alpha[i]];
        }
    }
}


trace_t
*make_empty_trace(void) {
     trace_t *R;
     R = (trace_t*)malloc(sizeof(*R));
     assert(R!=NULL);
     R->head = R->tail = NULL;
     return R;
 }

trace_t
*insert_at_tail(trace_t* R, action_t* addr) {
     step_t *new;
     new = (step_t*)malloc(sizeof(*new));
     assert(R!=NULL && new!=NULL);
     new->action = addr;
     new->next = NULL;
     if (R->tail==NULL) { /* this is the first insertion into the trace */
         R->head = R->tail = new;
     } else {
         R->tail->next = new;
         R->tail = new;
     }
     return R;
 }

void
free_trace(trace_t* R) {
     step_t *curr, *prev;
     assert(R!=NULL);
     curr = R->head;
     while (curr) {
         prev = curr;
         curr = curr->next;
         free(prev);
    }
    free(R);
 }

void 
print_out_state(state_t state) {
    int m;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    for (m=0;m<ASIZE;m++) {
        printf("%d", state[(int)alpha[m]]);
    }
    printf("\n");
}

void 
build_state(char* str, state_t state) {
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    int i;
    for (i=0;i<ASIZE;i++) {
        state[(int)alpha[i]] = 0;
    }
    assert(str!=NULL);
    /*update the initial states into the state array*/   
    for (i=0;str[i]!='\0';i++) {
        if (isalpha(str[i])) {
            state[(int)str[i]] = ON;
        }
    }
}

void 
turn_off(char* str, state_t* state) {
    int i;
    assert(str!=NULL);
    
    /*update the initial states into the state array*/   
    for (i=0;str[i]!='\0';i++) {
        if (isalpha(str[i])) {
            (*state)[(int)str[i]] = TURNOFF;
        }
    }
    // 2 means turn off this! 
}

/*this function reads action and stores each line into action_t*/
action_t* 
read_action(char str[], int start, int end) {
    action_t* new_action = (action_t*)malloc(sizeof(*new_action));
    char* pre;
    char* eff;
    char* turnoff;
    char* pre_off;
    int i, j=0, k, l=0, n=0, h, tf=0, preoff, po=0;
    
    pre = (char*)malloc(sizeof(*pre));
    assert(pre!=NULL);
    eff = (char*)malloc(sizeof(*eff));
    assert(eff != NULL);
    turnoff = (char*)malloc(sizeof(*turnoff));
    assert(turnoff != NULL);
    pre_off = (char*)malloc(sizeof(*pre_off));
    assert(pre_off != NULL);
    
    /*building preconditions*/
    for (i=start;str[i]!=':';i++) {
        pre[j++] = str[i];
    }
    pre[j] = '\0';
    
    for (preoff = i+1; str[preoff]!=':';preoff++) {
        pre_off[po++] = str[preoff];
    }
    pre_off[po] = '\0';
    
    // build precons then turn it off!
    build_state(pre, new_action->precon);
    turn_off(pre_off, &new_action->precon);
   
 
    
    for (i=start;i<end;i++) {
       if (isupper(str[i])) {
           new_action->name = str[i];
           n = i+2;
           break;
       }
    }
    
    
    for (k=n;str[k]!=':';k++) {
        eff[l++] = str[k];
    }
    eff[l] ='\0';
    build_state(eff, new_action->effect);
    
    for (h=k+1;h<end;h++) {
        turnoff[tf++] = str[h];
    }
    turnoff[tf] = '\0';
    turn_off(turnoff, &new_action->effect);
    /*a free for each malloc!!*/
    free(pre);
    pre=NULL;
    free(turnoff);
    turnoff = NULL;
    free(eff);
    eff=NULL;
    free(pre_off);
    pre_off = NULL;
    return new_action;
}

// this function accepts a precondition and a state, and see if the 'state' 
// satisfies the precondition
int 
check_precon(state_t pre, state_t state) {
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    int i;
    for (i=0;i<ASIZE;i++){
        if ((pre[(int)alpha[i]] == ON) && (state[(int)alpha[i]]== OFF)) {
            return 0;
        }
        else if ((pre[(int)alpha[i]] == ON) && 
                 (state[(int)alpha[i]]== TURNOFF)) {
            return 0;
        }
        else if ((pre[(int)alpha[i]] == TURNOFF) && 
                  (state[(int)alpha[i]]== ON)) {
            return 0;
        }
    }
    return 1;
}

/*this function implements effects of each unique command*/
void 
do_effect(state_t *effect, state_t *state) {
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    int i;
    for (i=0;i<ASIZE;i++) {
        if (((*effect)[(int)alpha[i]] == ON) && 
            ((*state)[(int)alpha[i]]== OFF)) {
            (*state)[(int)alpha[i]] = ON;
        }
        /*turn off this command if it is 2*/
        else if (((*effect)[(int)alpha[i]] == TURNOFF) && 
                ((*state)[(int)alpha[i]]== ON)) {
            (*state)[(int)alpha[i]] = OFF;
        }
    } 
    
} 
/*this function turnes on and off according to the 'effect' array received,
and modifies the 'state' array*/
void 
turn_it_off(state_t *effect, state_t *state) {
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    int i;
    for (i=0;i<ASIZE;i++) {
        if ((*effect)[(int)alpha[i]] == ON) {
            (*state)[(int)alpha[i]] = ON;
        }
        else if ((*effect)[(int)alpha[i]] == TURNOFF) {
            (*state)[(int)alpha[i]] = TURNOFF;
        }
    } 
    
} 

 
void 
exit_if_null(void *ptr, char *msg) {
	if (!ptr) {
		printf("unexpected null pointer: %s\n", msg);
		exit(EXIT_FAILURE);
	}
}

int 
is_same(state_t A, state_t B) {
    int m;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    for (m=0;m<ASIZE;m++) {
        if (A[(int)alpha[m]] != B[(int)alpha[m]]) return 0;
    }
    return 1;
}

void 
copy_states(state_t source, state_t dest) {
    int m;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    for (m=0;m<ASIZE;m++) {
        dest[(int)alpha[m]] = source[(int)alpha[m]];
    }
}

int 
is_modified(state_t effect, state_t state) {
    int m;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    for (m=0;m<ASIZE;m++) {
        if (effect[(int)alpha[m]] == ON && state[(int)alpha[m]] == 0) return 1;
    }
    return 0;
}
 
int 
search_for_stagetwo(trace_t* trace, state_t copy_state, char* commands, 
                    int cm, int start_pos, state_t variables_state) {

    step_t* tmp;
    state_t state, new, off;
    step_t* start_step;
    char* alpha = "abcdefghijklmnopqrstuvwxyz";
    int i=0, k=0, j;
    assert(trace->head!=NULL);
    /*make a copy state first*/
    for (i=0;i<ASIZE;i++) {
        state[(int)alpha[i]] = OFF;
    }
    for (i=0;i<ASIZE;i++) {
        off[(int)alpha[i]] = OFF;
    }

    // new is a state indicating variables that have been modified
    for (j=0;j<ASIZE;j++){
        if (variables_state[(int)alpha[j]] == TURNOFF){
            new[(int)alpha[j]] = 1;
        }
        else {
            new[(int)alpha[j]] = variables_state[(int)alpha[j]];
        }
    }
    tmp = trace->head;
    for (i=0;i<start_pos;i++){
        tmp = tmp->next;
    } 
  
    /*step through each command in trace, then do effect, check if same
    check if it satisfies the precon of the first element*/
    start_step = tmp;
    while (tmp != NULL) {  
        do_effect(&tmp->action->effect, &state);
        turn_it_off(&tmp->action->effect, &off);
        k++;
        
        if (is_same(state, copy_state)) {
            if (is_modified(tmp->action->effect, new) && 
                check_precon(start_step->action->precon, copy_state)) {
                // need a function to check if conditions are set back 
                printf("%5d: ", start_pos);
                print_int_array(commands, start_pos, start_pos+k); 
                return start_pos+k;
            }
            else {
                printf("%5d: ", start_pos);
                print_int_array(commands, start_pos, start_pos+k); 
                return start_pos+k;
            }
        }
        tmp = tmp->next; 
    }   
    return start_pos+1;
}

void 
print_int_array(char str[], int start, int end) {
    int i;
    for (i=start;i<end;i++) {
        printf("%c", str[i]);
    }
    printf("\n");
}

/*Algorithms are fun!!!*/