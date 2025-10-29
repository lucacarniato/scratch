import torch
import torch.nn.functional as F

def eikonal_grid_loss(sdf_grid: torch.Tensor, voxel_size, mask: torch.Tensor = None, eps: float = 1e-9):
    """
    sdf_grid:  (B,1,D,H,W) predicted SDF volume (PyTorch tensor)
    voxel_size: float or (vz, vy, vx) spacing in same units as SDF
    mask:      (B,1,D,H,W) optional weights (e.g., near-surface band)
    returns:   scalar loss = mean (||∇f|| - 1)^2 (masked mean if mask provided)
    """
    assert sdf_grid.ndim == 5 and sdf_grid.shape[1] == 1, "Expected shape (B,1,D,H,W)"
    if isinstance(voxel_size, (tuple, list)):
        vz, vy, vx = float(voxel_size[0]), float(voxel_size[1]), float(voxel_size[2])
    else:
        vz = vy = vx = float(voxel_size)

    # Central differences with replicate padding to preserve shape
    dz = (F.pad(sdf_grid[:,:,2:,:,:], (0,0,0,0,0,2), mode='replicate')
          - F.pad(sdf_grid[:,:,:-2,:,:], (0,0,0,0,2,0), mode='replicate')) / (2.0 * vz)
    dy = (F.pad(sdf_grid[:,:,:,2:,:], (0,0,0,2,0,0), mode='replicate')
          - F.pad(sdf_grid[:,:,:,:-2,:], (0,0,2,0,0,0), mode='replicate')) / (2.0 * vy)
    dx = (F.pad(sdf_grid[:,:,:,:,2:], (0,2,0,0,0,0), mode='replicate')
          - F.pad(sdf_grid[:,:,:,:,:-2], (2,0,0,0,0,0), mode='replicate')) / (2.0 * vx)

    grad_norm = torch.sqrt(dx*dx + dy*dy + dz*dz + eps)  # (B,1,D,H,W)
    eik = (grad_norm - 1.0) ** 2

    if mask is not None:
        m = mask.to(dtype=eik.dtype, device=eik.device)
        denom = torch.clamp(m.sum(), min=1.0)
        return (eik * m).sum() / denom
    else:
        return eik.mean()
