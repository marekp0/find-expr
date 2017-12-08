#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <mutex>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

constexpr u32 MAX_TERMS = 9;

// the types of operations we can do
enum Op : u32 {
    Concat,
    Add,
    Sub,
    Mul,
    Div,
    Exp,

    NUM_OPS
};

char opToChar(Op op)
{
    static char opMap[] = {
        '\0',
        '+',
        '-',
        '*',
        '/',
        '^'
    };
    return opMap[op];
}

// small helper to calculate x^y for integers
u32 ipow(u32 x, u32 y)
{
    u32 r = 1;
    for (u32 i = 0; i < y; i++)
        r *= x;

    return r;
}

class OpList {
public:
    std::vector<Op> ops;

    OpList(u32 size, u32 num)
    {
        ops.resize(size);

        // initialize, check bounds
        u32 d = ipow(NUM_OPS, size - 1);
        if (num >= d * NUM_OPS) {
            printf("OpList constructor: num too big!\n");
            exit(1);
        }

        // essentially convert num into a "size"-digit base-"num" integer
        for (u32 i = 0; i < size; i++) {
            ops[i] = (Op)(num / d);
            num = num % d;
            d /= NUM_OPS;
        }
    }

    // goes to the next possible set of operations
    void next()
    {
        for (u32 j = 0; j < ops.size(); j++) {
            // increment from right to left
            u32 i = ops.size() - j - 1;

            ops[i] = (Op)((u32)ops[i] + 1);
            if (ops[i] == NUM_OPS) {
                ops[i] = (Op)0;
            }
            else {
                break;
            }
        }
    }
};

// permutations of all possible orders of operations for various numbers of operations
// e.g. orderPerms[8] should contain all permutations of orders of 8 operations
//std::vector<std::vector<std::vector<u32>>> orderPerms;

typedef std::vector<u32> OpOrder;
typedef std::vector<OpOrder> OpOrderList;

void orderPermImpl(const OpOrder& order, u32 pos, OpOrderList& result)
{
    if (pos == order.size() - 1) {
        // end case: nothing more to generate, add it to the result
        result.push_back(order);
    }
    else {
        // swap the element as position pos with all others, the recursively
        // call this function
        for (u32 i = 0; i < order.size() - pos; i++) {
            std::vector<u32> newOrder = order;
            std::swap(newOrder[pos], newOrder[pos + i]);
            orderPermImpl(newOrder, pos + 1, result);
        }
    }
}

OpOrderList genOpOrders(u32 size)
{
    OpOrderList perms;
    OpOrder start;

    // create a vector of 0, 1, ..., size-1
    for (u32 i = 0; i < size; i++) {
        start.push_back(i);        
    }
    
    // generate all permutations of "start" and store them in "perms"
    orderPermImpl(start, 0, perms);

    return perms;   // move semantics hopefully?
}

struct Instruction {
    u32 inLeft, inRight;    // input "registers"
    u32 out;                // output "register"
    u32 opId;               // index into the list of operations
};

typedef std::vector<Instruction> Program;
typedef std::vector<Program> ProgramList;

// collection of programs for all numbers of operations, and all permutations of
// operation orders
std::vector<ProgramList> programs;

Program createProgram(const OpOrder& order)
{
    Program prog;

    // Keep track of register mappings. Initially, each register maps to itself.
    // As we do some operations, registers will map to other registers to
    // implement temporary values.
    std::vector<u32> regMap(order.size() + 1);
    for (int i = 0; i < regMap.size(); i++)
        regMap[i] = i;
    
    for (auto i : order) {
        Instruction insn;
        insn.inLeft = regMap[i];
        insn.inRight = regMap[i + 1];
        // store result in the lower of the two input registers
        insn.out = std::min(insn.inLeft, insn.inRight);
        insn.opId = i;

        prog.push_back(insn);

        // map all of the input registers to the output register
        for (auto& m : regMap) {
            if (m == insn.inLeft || m == insn.inRight)
                m = insn.out;
        }
    }

    // the final result should be in register 0 at the end

    return prog;
}

// list of all programs with 8 operations
ProgramList allProgs;

// generates allProgs
void genAllProgs(u32 numOps)
{
    auto opOrders = genOpOrders(numOps);
    for (auto& o : opOrders) {
        allProgs.push_back(createProgram(o));
    }
}

// debug
void dumpProg(const Program& p)
{
    for (auto i : p) {
        printf("%d: %d, %d --> %d\n", i.opId, i.inLeft, i.inRight, i.out);
    }
}

struct Expression {
    std::vector<double> terms;
    OpList ol;
    Program* prog;

    bool isValid()
    {
        for (u32 i = 0; i < prog->size(); i++) {
            Instruction& insn = (*prog)[i];
            Op op = ol.ops[insn.opId];
            if (op == Concat && i > 0) {
                Instruction& prevInsn = (*prog)[i - 1];
                // ensure that all concats are at the beginning
                if (ol.ops[prevInsn.opId] != Concat)
                    return false;

                // ensure that the order is correct
                if (insn.inRight <= prevInsn.inLeft)
                    return false;
            }
        }
        return true;
    }

    // evaluates the expression
    double eval()
    {
        double reg[MAX_TERMS]; 

        // initialize registers
        for (int i = 0; i < terms.size(); i++) {
            reg[i] = terms[i];
        }

        for (auto& insn : *prog) {
            // get input values and operation
            double left = reg[insn.inLeft];
            double right = reg[insn.inRight];
            Op op = ol.ops[insn.opId];

            // calculate output (R.I.P. branch predictor)
            double out;
            switch (op) {
            case Concat:
                // assume that all concatenations are done at the beginning, and
                // in order from left to right
                out = left * 10.0 + right;
                break;
            case Add:
                out = left + right;
                break;
            case Sub:
                out = left - right;
                break;
            case Mul:
                out = left * right;
                break;
            case Div:
                out = left / right;
                break;
            case Exp:
                out = pow(left, right);
                break;
            default:
                printf("We weren't supposed to get here!\n");
                exit(1);
            }

            // store output
            reg[insn.out] = out;
        }

        // result should be in the first register
        return reg[0];
    }

    void dump()
    {
        for (u32 i = 0; i < terms.size(); i++) {
            printf("%d", (int)terms[i]);
            if (i != terms.size() - 1)
                printf("%c", opToChar(ol.ops[i]));
        }
        printf("; ");
        for (auto& insn : *prog) {
            printf("%d ", insn.opId);
        }
        printf("\n");
    }
};

std::mutex mtxWrite;

void tryExpressions(std::vector<double> terms, u32 start, u32 end, double goal)
{
    u32 numIterations = end - start;

    Expression e {
        terms,
        OpList(terms.size() - 1, start),
        nullptr
    };

    for (u32 i = 0; i < numIterations; i++) {
        for (auto& p : allProgs) {
            e.prog = &p;
            if (e.isValid()) {
                double result = e.eval();
                if (fabs(result - goal) < 0.01) {
                    // lock the mutex to make sure multiple lines from other
                    // threads aren't clobbered together
                    mtxWrite.lock();
                    e.dump();
                    mtxWrite.unlock();
                }
            }
        }

        e.ol.next();   
    }
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        printf("Usage: findExpr <num>");
        return 1;
    }
    u32 GOAL = atoi(argv[1]);

    // generate all programs
    genAllProgs(MAX_TERMS - 1);

    // generate terms
    std::vector<double> terms;
    for (u32 i = 0; i < MAX_TERMS; i++)
        terms.push_back(i + 1);

    // create threads
    u32 numThreads = std::thread::hardware_concurrency();
    u32 numOpLists = ipow(NUM_OPS, MAX_TERMS - 1);
    u32 stride = numOpLists / numThreads;
    std::vector<std::thread*> threads;
    for (u32 i = 0; i < numThreads; i++) {
        threads.push_back(new std::thread(tryExpressions, terms, i * stride, (i + 1) * stride, GOAL));
    }

    for (auto t : threads) {
        t->join();
        delete t;
    }
    printf("done!\n");
    return 0;
}

