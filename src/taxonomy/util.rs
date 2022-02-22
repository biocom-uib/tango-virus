use super::NodeId;

#[derive(Clone)]
pub struct ParentsIter<GetParent> {
    get_parent: GetParent,
    current: NodeId,
}

impl<GetParent> ParentsIter<GetParent> {
    pub fn new_with(node_id: NodeId, get_parent: GetParent) -> Self {
        Self {
            current: node_id,
            get_parent: get_parent,
        }
    }
}

impl<GetParent> Iterator for ParentsIter<GetParent>
where
    GetParent: FnMut(NodeId) -> Option<NodeId>,
{
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(parent) = (self.get_parent)(self.current) {
            self.current = parent;
            Some(parent)
        } else {
            // root
            None
        }
    }
}


#[derive(Clone)]
pub struct PreOrderIter<GetChildren, ChildrenIter> {
    get_children: GetChildren,
    stack: Vec<ChildrenIter>,
}

impl<GetChildren, ChildrenIter> PreOrderIter<GetChildren, ChildrenIter>
where
    GetChildren: FnMut(NodeId) -> ChildrenIter,
{
    pub fn new_with(node_id: NodeId, mut get_children: GetChildren) -> Self {
        Self {
            stack: vec![get_children(node_id)],
            get_children: get_children,
        }
    }
}

impl<GetChildren, ChildrenIter> Iterator for PreOrderIter<GetChildren, ChildrenIter>
where
    GetChildren: FnMut(NodeId) -> ChildrenIter,
    ChildrenIter: Iterator<Item = NodeId>,
{
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(node) = self.stack.last_mut() {
            if let Some(child) = node.next() {
                self.stack.push((self.get_children)(child));
                Some(child)
            } else {
                self.stack.pop();
                self.next()
            }
        } else {
            None
        }
    }
}


#[derive(Clone)]
pub struct PostOrderIter<GetChildren, ChildrenIter: Iterator> {
    get_children: GetChildren,
    stack: Vec<(NodeId, ChildrenIter)>,
}

impl<GetChildren, ChildrenIter> PostOrderIter<GetChildren, ChildrenIter>
where
    GetChildren: FnMut(NodeId) -> ChildrenIter,
    ChildrenIter: Iterator<Item = NodeId>,
{
    pub fn new_with(node_id: NodeId, mut get_children: GetChildren) -> Self {
        Self {
            stack: vec![(node_id, get_children(node_id))],
            get_children: get_children,
        }
    }
}

impl<GetChildren, ChildrenIter> Iterator for PostOrderIter<GetChildren, ChildrenIter>
where
    GetChildren: FnMut(NodeId) -> ChildrenIter,
    ChildrenIter: Iterator<Item = NodeId>,
{
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(child) = children.next() {
                let grand_children = (self.get_children)(child);
                self.stack.push((child, grand_children));
                self.next()
            } else {
                let current = *current;
                self.stack.pop();
                Some(current)
            }
        } else {
            None
        }
    }
}
