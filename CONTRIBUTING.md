# Contributing to Mgart

Contributions are very welcome! 
Please feel free to contribute anything you got, be that simple 
questions, a friendly chat about algorithmic art, feature requests, 
bug reports, bug fixes, or even implementations of new algorithms. 

  
## Development Guidelines

There are many great tools to generate digital art with.
No software has been able to deliver an all encumbering experience
when it comes to the generation and exploration of algorithmic 
beauty, algorithm and hardware support, or speed.
Mgart definetly does not aim to be the tool that changes that.
The goal of Mgart is more modest, to be a tool to generate beautiful 
digital imagery with, using a simple, declarative 
[DSL](https://en.wikipedia.org/wiki/Domain-specific_language).
When you want to contribue source code, please follow these 
development guidelines:

* **API-first:** The CLI and DSL should be simple, effective and 
  self-explanatory.

* **Code quality over performance:** Quality code that is well
  abstracted, well tested and readable is valued more than code
  that is slightly faster but unreadable.

* **Correctness violations:** Mgart tries to, but does not 
  necessarily produce images containing mathematically correct  
  renderings of an algorithm. 
  Violations of the mathematical constraints in favor of speed or 
  beauty must be specified by the user in the DSL.

* **Fast algorithms over fast hardware:** Rather than adding support
  for various GPUs that perform brute-force faster than a CPU, Mgart
  tries to implement smart algorithms that reduce the amount of 
  computation.


## Code of Conduct

Mgart adheres to the [Contributor Covenant v2.1](CODE_OF_CONDUCT.md).
