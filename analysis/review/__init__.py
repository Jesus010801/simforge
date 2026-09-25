"""Review-session preparation (Trajectory Review Mode, headless).

Sessionizes the scientific backend in ``analysis.campaign`` into one
reproducible, reopenable :class:`~analysis.review.dataset.ReviewDataset`.
No scientific definitions live here, and nothing starts a viewer or server.
"""
from analysis.review.dataset import (  # noqa: F401
    DATASET_FILENAME, REVIEW_SCHEMA_VERSION, Availability, CapabilityState, ObservableEntry,
    ReviewDataset, ReviewResult,
)
from analysis.review.prepare import Preparation, ReviewError, prepare_review  # noqa: F401
from analysis.review.request import (  # noqa: F401
    DisplayRequest, ObservableRequest, RequestError, ReviewRequest, load_request_file,
    parse_show,
)
from analysis.review.validate import (  # noqa: F401
    ReviewValidation, load_review_dataset, open_review_session, validate_review_dataset,
)
