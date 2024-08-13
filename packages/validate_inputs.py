def validate_inputs(x, y, z)-> tuple[str, str, str]:
    # Convert each input to lowercase
    x = x.lower()
    y = y.lower()
    z = z.lower()
        
    # Ensure each variable is of length one
    assert len(x) == 1, "x must be of length 1"
    assert len(y) == 1, "y must be of length 1"
    assert len(z) == 1, "z must be of length 1"
    
    # Ensure each variable is one of 'a', 'b', or 'c'
    allowed_values = {"x", "y", "z"}
    assert x in allowed_values, "x must be one of 'x', 'y', or 'z'"
    assert y in allowed_values, "y must be one of 'x', 'y', or 'z'"
    assert z in allowed_values, "z must be one of 'x', 'y', or 'z'"
    
    # Ensure all variables are unique
    assert len({x, y, z}) == 3, "x, y, and z must be unique"

    print("All inputs are valid.")
    
    return x, y, z

# # Example usage:
# validate_inputs("x", "y", "z")  # This will pass validation.
# validate_inputs("x", "y", "y")  # This will raise an AssertionError.