Sphinx is used for the documentation of this package.
However, Sphinx does not know about Abaqus package and ignores modules,
if a package cannot be imported.
To prevent this, `doc/source/dummy_packages` contains dummy python files that look like Abaqus packages.
Sphinx can import these dummy python files and hence generate a documentation.

Follow these rules to keep the Sphinx documentation running:

- import Abaqus modules using the syntax:
  ```python
  import abaqus
  import abaqusConstants
  import odbAccess
  ```
  Aliases may be used such as `import abaqus as abq`.
- **Do not import specific functions**, as they are not defined in the dummy python files:
  ```python
  from abaqus import mdb
  ```
- Call Abaqus specific code only in functions or in the main-clause, but not in the global scope.
  - **NOT OK**
    ```python
    import abaqus as abq
    
    abq.Mdb()  # this is the global scope that is executed by importing a module
    ```
  - **OK**
    ```python
    import abaqus as abq
    
    if __name__ == '__main__':
        abq.Mdb()  # main clause is not executed during the import
    ```
  - **OK**
    ```python
    import abaqus as abq
    
    def some_computation():
        abq.Mdb()  # not executed, as long as the function is not called in the outer scope.
    ```