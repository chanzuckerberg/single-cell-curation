# Current state
The current validator uses a yaml file which defines the schema and a python validator that interprets it. Over time, the validation requirements have grown in complexity. Unfortunately the schema defined in the yaml has not evolved to keep up. This has resulted in complexity being added to the validator when it should be in the yaml. This has made the validator difficult to test and maintain.


There is a severe lack of unit tests for the validator. This is because the validator is difficult to test. The specific rules to test are not well isolated and buried in if-else logic trees.
The current setup is good, but is not well suited for testing and it is reaching its limit in complexity.

# Yaml
- element_type key word uses parsing logic hidden in the validator. It could be written more clearly and reusably.
- complex_rule and it's keys should be represented as a function with its keys as parameters. This would help make testing it easier.
- deprecated_columns and reserved_columns should not be used in type: dict
- we have a schema that validates and transforms the data which means its user for both validating the anndata and writing the labels.
```    
index:
  unique: true
  type: feature_id
  add_labels:
    - type: feature_id
```
- The yaml schema should be highly readable and configurable, where as the validation code should be written as performant and testable.

## Problem with new tools
they will all require thinking about how we will translate the transformation keywords in the yaml. The libraries do support modifying the model at runtime, which means we won't be able to run validation and adding labels separately. More investigation is needed to see how we can work around this.

# Cerberus
- License: ISC
- Website: https://docs.python-cerberus.org/en/stable/
- 3.1K stars on github


Cerberus is a lightweight and extensible data validation library for Python. Cerberus provides type checking and other validation features to ensure that the data is in the expected format.


## Pros
- will be able to remove a lot of the if logic in the validator. This will make the code easier to read and test.
- Gives us access to [validation workflow](https://docs.python-cerberus.org/validation-rules.html) features we would have to implement ourselves.
- unit testing would be limited to verifying customer rules and checks are working
- supports checks that depend on other fields using [Validator.document](https://docs.python-cerberus.org/customize.html#validator-document)
- allows us to continue to use the yaml schema and expand on it.
- support custom [rules and checks](https://docs.python-cerberus.org/customize.html#custom-rules)
- many if the validation rules already follow the cerberus naming convention for example `_check_deprecated_columns()`

## Cons
- will need to rewrite the yaml schema to work with cerberus
- will need to break down the anndata object into its components to validate each component because Cerberus assumes everything is a Mapping.
- We will need to add support for dataframes and numpy arrays and any other nonstandard data types.
- only throws errors, no off the shelf support for warnings. We'd need to subclass the error handler to support warnings.


# pydantic
license: MIT
website: https://pydantic-docs.helpmanual.io/


lots of community support
18.5k stars on github
lots of examples


## Pros
- very large community
- works with pandera for dataframes
- it supports parallelized validation of models.
- support validating class objects directly. Don't need to break down the anndata object into its components, but breaking it down will make it easier to test.
- support for custom [validators](https://pydantic-docs.helpmanual.io/usage/validators/)


## Cons
- would need a full rewrite of the validator
- would no longer use the yaml schema
- ease of testing would depend on how well the anndata object can be broken down into its components
- Passing context between models is not straightforward which makes rules that depend on other fields difficult to implement.
- not as straightforward for a non python developer to understand the validation rules.
- steeper learning curve for the team.


# pandera
license: MIT
website: https://www.union.ai/pandera
3K stars on github


integrates with pydantic, mainly used for dataframes.


If pydantic is used then we mine as well use pandera to validate the dataframes.


# great expectations
license: Apache 2.0
website: https://greatexpectations.io/
9.4K stars on github


This is a large chain of tools for data validation. I did not have time to evaluate it, but based on the internet it would not be trivial to migrate to this tool, but it would be a good long term solution.


# Links
- [comparing pandera and great expectations](https://endjin.com/blog/2023/03/a-look-into-pandera-and-great-expectations-for-data-validation)

# Next steps
What ever course of action we take, we should start by writing a test suite for the current validator. This will help us understand the current behavior and make sure we don't break anything when we migrate to a new tool.


# Optimzing the current validator
- only check for canonical sparse array during seurat conversion
- remove the need for the cxg converter to check for the sparsity of the matrix by returning the sparsity from the validator
  - Use a single sparsity check for both the validator and the cxg converter
