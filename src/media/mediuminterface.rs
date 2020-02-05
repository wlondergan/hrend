/// Distributing light through a medium requires knowledge of the inside and outside of the medium.
/// The MediumInterface struct holds information about the media contained within it to make this
/// interaction easier to manage.
#[derive(Copy, Clone, PartialEq)]
pub struct MediumInterface { //TODO implement
   // inside: Medium,
   // outside: Medium
}

impl MediumInterface {
    pub fn default() -> MediumInterface {
        MediumInterface::new()
    }

    pub fn new() -> MediumInterface {
        MediumInterface {
            
        }
    }
}