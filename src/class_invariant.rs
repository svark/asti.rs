pub trait ClassInvariant
{
    fn is_valid(&self) -> Result<bool, &str>;
}
