use super::*;

#[test]
fn test_stays_inline_up_to_capacity() {
    let list: InlineOrHeapList<u32, 3> = [1, 2, 3].into_iter().collect();
    assert!(matches!(list, InlineOrHeapList::Inline(_)));
    assert_eq!(*list, [1, 2, 3]);
}

#[test]
fn test_moves_to_the_heap_past_capacity_and_keeps_order() {
    let list: InlineOrHeapList<u32, 2> = [1, 2, 3, 4, 5].into_iter().collect();
    assert!(matches!(list, InlineOrHeapList::Heap(_)));
    assert_eq!(*list, [1, 2, 3, 4, 5]);
}

#[test]
fn test_equality_ignores_the_storage() {
    let inline: InlineOrHeapList<u32, 4> = [3, 1].into_iter().collect();
    let heap: InlineOrHeapList<u32, 4> = vec![3, 1].into();
    assert_eq!(inline, heap);
    assert_eq!(format!("{inline:?}"), format!("{:?}", vec![3, 1]));
}

#[test]
fn test_sorts_in_place_through_the_slice() {
    let mut list: InlineOrHeapList<u32, 4> = [3, 1, 2].into_iter().collect();
    list.sort_unstable();
    assert_eq!(*list, [1, 2, 3]);
    assert!(InlineOrHeapList::<u32, 4>::new().is_empty());
}
