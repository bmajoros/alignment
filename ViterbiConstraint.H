#ifndef INCL_ViterbiConstraint_H
#define INCL_ViterbiConstraint_H


enum ViterbiConstraint {
  VC_UNCONSTRAINED, // choose indel history and parses for parent & child
  VC_PARENT_PARSE,  // choose everything except parent parse
  VC_CHILD_PARSE,   // choose everything except child parse
  VC_INDEL_HISTORY, // choose only parent & child parses for fixed indel hist
  VC_INDEL_AND_PARENT_PARSE,// may choose only a parse for the child
  VC_INDEL_AND_CHILD_PARSE, // may choose only a parse for the parent
  VC_PARENT_AND_CHILD_PARSE // may choose an indel history; parses are fixed
};


#endif
