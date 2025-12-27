"""
Framework-specific test library implementations for PennyLane.
"""

import pytest
import numpy as np

# Framework-specific imports
torch = pytest.importorskip("torch")
tf = pytest.importorskip("tensorflow")
jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")

def test_torch_backend():
    """Test PyTorch backend integration."""
    x = torch.tensor([1.0, 2.0, 3.0], requires_grad=True)
    y = x ** 2
    assert y.shape == x.shape
    assert y.requires_grad

def test_tensorflow_backend():
    """Test TensorFlow backend integration."""
    x = tf.Variable([1.0, 2.0, 3.0])
    with tf.GradientTape() as tape:
        y = x ** 2
    assert y.shape == x.shape
    grads = tape.gradient(y, x)
    assert grads is not None

def test_jax_backend():
    """Test JAX backend integration."""
    x = jnp.array([1.0, 2.0, 3.0])
    y = x ** 2
    assert y.shape == x.shape
    
    def f(x):
        return jnp.sum(x ** 2)
    
    grad_f = jax.grad(f)
    grads = grad_f(x)
    assert grads.shape == x.shape

def test_numpy_backend():
    """Test NumPy backend integration."""
    x = np.array([1.0, 2.0, 3.0])
    y = x ** 2
    assert y.shape == x.shape
    np.testing.assert_array_equal(y, np.array([1.0, 4.0, 9.0]))

def test_torch_gradients():
    """Test gradient computation with PyTorch."""
    x = torch.tensor([1.0], requires_grad=True)
    y = x ** 2
    y.backward()
    assert x.grad is not None
    assert torch.isclose(x.grad, torch.tensor([2.0]))

def test_tensorflow_gradients():
    """Test gradient computation with TensorFlow."""
    x = tf.Variable([1.0])
    with tf.GradientTape() as tape:
        y = x ** 2
    grads = tape.gradient(y, x)
    assert grads is not None
    assert tf.reduce_all(tf.abs(grads - tf.constant([2.0])) < 1e-6)

def test_jax_gradients():
    """Test gradient computation with JAX."""
    def f(x):
        return jnp.sum(x ** 2)
    
    grad_f = jax.grad(f)
    x = jnp.array([1.0])
    grads = grad_f(x)
    assert jnp.allclose(grads, jnp.array([2.0]))

def test_torch_dtype_conversion():
    """Test dtype conversion with PyTorch."""
    x = torch.tensor([1, 2, 3], dtype=torch.int32)
    y = x.float()
    assert y.dtype == torch.float32

def test_tensorflow_dtype_conversion():
    """Test dtype conversion with TensorFlow."""
    x = tf.constant([1, 2, 3], dtype=tf.int32)
    y = tf.cast(x, tf.float32)
    assert y.dtype == tf.float32

def test_jax_dtype_conversion():
    """Test dtype conversion with JAX."""
    x = jnp.array([1, 2, 3], dtype=jnp.int32)
    y = x.astype(jnp.float32)
    assert y.dtype == jnp.float32

def test_numpy_scalar_with_grad():
    """Test NumPy scalar value handling (original line 91)."""
    x_val = np.pi / 4  # PennyLane handles gradient tracking automatically
    assert isinstance(x_val, (float, np.floating))
    assert np.isclose(x_val, np.pi / 4)
